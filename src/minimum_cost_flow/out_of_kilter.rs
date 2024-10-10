use crate::minimum_cost_flow::csr::CSR;
use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::status::Status;
use num_traits::NumAssign;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Neg;

// O(nU * (m + n) log n)
#[derive(Default)]
pub struct OutOfKilter<Flow> {
    csr: CSR<Flow>,
}

impl<Flow> OutOfKilter<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn solve(&mut self, graph: &mut Graph<Flow>) -> Status {
        if graph.is_unbalance() {
            return Status::Unbalanced;
        }

        let (_source, artificial_nodes, artificial_edges) = graph.construct_extend_network_feasible_solution();
        self.csr.build(graph);

        let mut out_of_kilter_edges = Vec::new();
        for (edge_id, edge) in self.csr.inside_edge_list.iter().enumerate() {
            let p = self.csr.inside_edge_list[edge.rev].to;
            if self.kilter_number(p, edge_id) != Flow::zero() {
                let q = edge.to;
                out_of_kilter_edges.push((p, q, edge_id));
            }
        }

        'outer: for (p, q, edge_id) in out_of_kilter_edges {
            while self.kilter_number(p, edge_id) > Flow::zero() {
                let (dist, prev) = self.shortest_path(q);
                if prev[p].is_none() {
                    break 'outer;
                }

                // update potentials
                for u in 0..self.csr.num_nodes {
                    if let Some(d) = dist[u] {
                        self.csr.potentials[u] -= d;
                    }
                }

                // update flow
                let edge = &self.csr.inside_edge_list[edge_id];
                if self.csr.reduced_cost(p, edge) < Flow::zero() {
                    self.update_flow_in_cycle(q, edge_id, prev);
                }
            }
        }

        self.csr.set_flow(graph);

        let status = if artificial_edges.iter().all(|&edge_id| graph.edges[edge_id].flow == Flow::zero()) {
            Status::Optimal
        } else {
            Status::Infeasible
        };
        graph.remove_artificial_sub_graph(&artificial_nodes, &artificial_edges);

        status
    }

    fn kilter_number(&self, u: usize, edge_id: usize) -> Flow {
        let edge = &self.csr.inside_edge_list[edge_id];
        if self.csr.reduced_cost(u, edge) >= Flow::zero() {
            Flow::zero()
        } else {
            edge.residual_capacity()
        }
    }

    fn shortest_path(&mut self, s: usize) -> (Vec<Option<Flow>>, Vec<Option<usize>>) {
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

            for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
                let edge = &self.csr.inside_edge_list[edge_id];
                if edge.residual_capacity() <= Flow::zero() {
                    continue;
                }

                let new_dist = d.0 + self.csr.reduced_cost(u, edge).max(Flow::zero());
                if dist[edge.to].is_none() || dist[edge.to].unwrap() > new_dist {
                    dist[edge.to] = Some(new_dist);
                    prev[edge.to] = Some(edge_id);
                    bh.push((Reverse(new_dist), edge.to));
                }
            }
        }

        (dist, prev)
    }

    fn update_flow_in_cycle(&mut self, q: usize, edge_id: usize, mut prev: Vec<Option<usize>>) {
        prev[q] = Some(edge_id); // p -> q

        // calculate delta
        let mut delta = self.csr.inside_edge_list[edge_id].residual_capacity();
        let mut v = q;
        while let Some(edge_idx) = prev[v] {
            delta = delta.min(self.csr.inside_edge_list[edge_idx].residual_capacity());
            let rev = self.csr.inside_edge_list[edge_idx].rev;
            v = self.csr.inside_edge_list[rev].to;
            if v == q {
                break;
            }
        }

        // update flow
        let mut v = q;
        while let Some(edge_id) = prev[v] {
            let rev = self.csr.inside_edge_list[edge_id].rev;
            v = self.csr.inside_edge_list[rev].to;
            self.csr.push_flow(v, edge_id, delta);
            if v == q {
                break;
            }
        }
    }
}
