use crate::minimum_cost_flow::csr::CSR;
use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::status::Status;
use num_traits::NumAssign;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Neg;

#[derive(Default)]
pub struct SuccessiveShortestPath<Flow> {
    csr: CSR<Flow>,
}

impl<Flow> SuccessiveShortestPath<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn solve(&mut self, graph: &mut Graph<Flow>) -> Status {
        if graph.is_unbalance() {
            return Status::Unbalanced;
        }
        self.csr.build(graph);

        for s in 0..self.csr.num_nodes {
            while self.csr.excesses[s] > Flow::zero() {
                match self.calculate_distance(s) {
                    Some((t, visited, dist, prev)) => {
                        // update potentials
                        for u in 0..self.csr.num_nodes {
                            if visited[u] {
                                self.csr.potentials[u] = self.csr.potentials[u] - dist[u].unwrap() + dist[t].unwrap();
                            }
                        }
                        // update flow
                        self.update_flow(s, t, prev);
                    }
                    None => break,
                }
            }
        }

        self.csr.set_flow(graph);

        if self.csr.excesses.iter().all(|&e| e == Flow::zero()) {
            Status::Optimal
        } else {
            Status::Infeasible
        }
    }

    pub fn calculate_distance(&mut self, s: usize) -> Option<(usize, Vec<bool>, Vec<Option<Flow>>, Vec<Option<usize>>)> {
        let mut prev = vec![None; self.csr.num_nodes];
        let mut bh = BinaryHeap::new();
        let mut dist: Vec<Option<Flow>> = vec![None; self.csr.num_nodes];
        let mut visited = vec![false; self.csr.num_nodes];

        bh.push((Reverse(Flow::zero()), s));
        dist[s] = Some(Flow::zero());

        while let Some((d, u)) = bh.pop() {
            if visited[u] {
                continue;
            }
            visited[u] = true;

            if self.csr.excesses[u] < Flow::zero() {
                return Some((u, visited, dist, prev));
            }

            for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
                let edge = &self.csr.inside_edge_list[edge_id];
                if edge.residual_capacity() == Flow::zero() {
                    continue;
                }

                let new_dist = d.0 + self.csr.reduced_cost(u, edge);
                if dist[edge.to].is_none() || dist[edge.to].unwrap() > new_dist {
                    dist[edge.to] = Some(new_dist);
                    prev[edge.to] = Some(edge_id);
                    bh.push((Reverse(new_dist), edge.to));
                }
            }
        }

        None
    }

    fn update_flow(&mut self, s: usize, t: usize, prev: Vec<Option<usize>>) {
        debug_assert!(self.csr.excesses[s] > Flow::zero() && self.csr.excesses[t] < Flow::zero());

        // calculate delta
        let mut delta = self.csr.excesses[s].min(-self.csr.excesses[t]);
        {
            let mut v = t;
            while let Some(edge_idx) = prev[v] {
                delta = delta.min(self.csr.inside_edge_list[edge_idx].residual_capacity());
                let rev = self.csr.inside_edge_list[edge_idx].rev;
                v = self.csr.inside_edge_list[rev].to;
            }
            delta = delta.min(self.csr.excesses[v]);
            debug_assert_eq!(s, v);
            debug_assert!(delta > Flow::zero());
        }

        // update flow
        {
            let mut v = t;
            while let Some(edge_idx) = prev[v] {
                // push
                let rev = self.csr.inside_edge_list[edge_idx].rev;
                self.csr.inside_edge_list[edge_idx].flow += delta;
                self.csr.inside_edge_list[rev].flow -= delta;
                v = self.csr.inside_edge_list[rev].to;
            }
            debug_assert_eq!(s, v);
        }

        self.csr.excesses[t] += delta;
        self.csr.excesses[s] -= delta;
    }
}
