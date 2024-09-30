use crate::minimum_cost_flow::csr::CSR;
use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::status::Status;
use num_traits::NumAssign;
use std::collections::{BinaryHeap, VecDeque};
use std::ops::Neg;

#[derive(Default)]
pub struct PrimalDual<Flow> {
    csr: CSR<Flow>,

    // maximum flow(dinic)
    que: VecDeque<usize>,
    distances: Vec<usize>,
    current_edge: Vec<usize>,
}

impl<Flow> PrimalDual<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn solve(&mut self, graph: &mut Graph<Flow>) -> Status {
        if graph.is_unbalance() {
            return Status::Unbalanced;
        }

        // transforms the minimum cost flow problem into a problem with a single excess node and a single deficit node.
        let (source, sink, artificial_nodes, artificial_edges) = graph.construct_extend_network_one_supply_one_demand();
        self.csr.build(graph);

        self.distances.resize(self.csr.num_nodes, 0);
        self.current_edge.resize(self.csr.num_nodes, 0);

        while self.csr.excesses[source] > Flow::zero() {
            if !self.dual(source, sink) {
                break;
            }
            self.primal(source, sink);
        }

        self.csr.set_flow(graph);

        graph.remove_artificial_sub_graph(&artificial_nodes, &artificial_edges);
        if self.csr.excesses[source] != Flow::zero() || self.csr.excesses[sink] != Flow::zero() {
            return Status::Infeasible;
        }

        Status::Optimal
    }

    // update potentials
    fn dual(&mut self, source: usize, sink: usize) -> bool {
        assert!(self.csr.excesses[source] > Flow::zero());

        // calculate the shortest path
        let mut dist: Vec<Option<Flow>> = vec![None; self.csr.num_nodes];
        let mut visited = vec![false; self.csr.num_nodes];
        {
            let mut bh: BinaryHeap<(Flow, usize)> = BinaryHeap::new();

            bh.push((Flow::zero(), source));
            dist[source] = Some(Flow::zero());

            while let Some((mut d, u)) = bh.pop() {
                d = -d;

                if visited[u] {
                    continue;
                }
                visited[u] = true;

                for edge_index in self.csr.start[u]..self.csr.start[u + 1] {
                    let e = &self.csr.inside_edge_list[edge_index];
                    if e.residual_capacity() == Flow::zero() {
                        continue;
                    }
                    if dist[e.to].is_none() || dist[e.to].unwrap() > d + self.csr.reduced_cost(u, e) {
                        dist[e.to] = Some(d + self.csr.reduced_cost(u, e));
                        bh.push((-dist[e.to].unwrap(), e.to));
                    }
                }
            }
        }

        // update potentials
        for u in 0..self.csr.num_nodes {
            if visited[u] {
                self.csr.potentials[u] -= dist[u].unwrap();
            }
        }

        visited[sink]
    }

    fn primal(&mut self, source: usize, sink: usize) {
        assert!(self.csr.excesses[source] > Flow::zero() && self.csr.excesses[sink] < Flow::zero());

        let mut flow = Flow::zero();
        while self.csr.excesses[source] > Flow::zero() {
            self.update_distances(source, sink);

            // no s-t path
            if self.distances[source] >= self.csr.num_nodes {
                break;
            }

            self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
            match self.dfs(source, sink, self.csr.excesses[source]) {
                Some(delta) => flow += delta,
                None => break,
            }
        }
        self.csr.excesses[source] -= flow;
        self.csr.excesses[sink] += flow;
    }

    // O(n + m)
    // calculate the distance from u to sink in the residual network
    // if such a path does not exist, distance[u] becomes self.num_nodes
    pub fn update_distances(&mut self, source: usize, sink: usize) {
        self.que.clear();
        self.que.push_back(sink);
        self.distances.fill(self.csr.num_nodes);
        self.distances[sink] = 0;

        while let Some(v) = self.que.pop_front() {
            for edge in self.csr.inside_edge_list[self.csr.start[v]..self.csr.start[v + 1]].iter() {
                // e.to -> v
                if edge.flow > Flow::zero() && self.distances[edge.to] == self.csr.num_nodes && self.csr.reduced_cost_rev(v, edge) == Flow::zero() {
                    self.distances[edge.to] = self.distances[v] + 1;
                    if edge.to != source {
                        self.que.push_back(edge.to);
                    }
                }
            }
        }
    }

    fn dfs(&mut self, u: usize, sink: usize, upper: Flow) -> Option<Flow> {
        if u == sink {
            return Some(upper);
        }

        let mut res = Flow::zero();
        for edge_index in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = edge_index;

            if !self.is_admissible_edge(u, edge_index) || self.csr.reduced_cost(u, &self.csr.inside_edge_list[edge_index]) != Flow::zero() {
                continue;
            }

            let v = self.csr.inside_edge_list[edge_index].to;
            let residual_capacity = self.csr.inside_edge_list[edge_index].residual_capacity();
            if let Some(d) = self.dfs(v, sink, residual_capacity.min(upper - res)) {
                let rev = self.csr.inside_edge_list[edge_index].rev;

                // update flow
                self.csr.inside_edge_list[edge_index].flow += d;
                self.csr.inside_edge_list[rev].flow -= d;

                res += d;
                if res == upper {
                    return Some(res);
                }
            }
        }
        self.current_edge[u] = self.csr.start[u + 1];
        self.distances[u] = self.csr.num_nodes;

        Some(res)
    }

    #[inline]
    pub fn is_admissible_edge(&self, from: usize, i: usize) -> bool {
        self.csr.inside_edge_list[i].residual_capacity() > Flow::zero() && self.distances[from] == self.distances[self.csr.inside_edge_list[i].to] + 1
    }
}
