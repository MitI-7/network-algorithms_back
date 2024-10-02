use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;
use std::collections::VecDeque;

#[derive(Default)]
pub struct CapacityScaling<Flow> {
    csr: CSR<Flow>,
    current_edge: Vec<usize>,
    que: VecDeque<usize>,
}

impl<Flow> CapacityScaling<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        self.csr.build(graph);
        self.current_edge.resize(self.csr.num_nodes, 0);
        let two = Flow::one() + Flow::one();

        let max_capacity = self.csr.inside_edge_list.iter().map(|e| e.upper).max().unwrap();
        let mut delta = Flow::one();
        while delta <= max_capacity {
            delta *= two;
        }
        delta /= two;

        let upper = self.csr.neighbors(source).fold(Flow::zero(), |sum, e| sum + e.upper);
        let mut flow = Flow::zero();
        while delta > Flow::zero() {
            // solve maximum flow in lambda-residual network
            loop {
                self.bfs(source, sink, delta);

                // no s-t path
                if self.csr.distances[source] >= self.csr.num_nodes {
                    break;
                }

                self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
                match self.dfs(source, sink, upper, delta) {
                    Some(delta) => flow += delta,
                    None => break,
                }
            }
            delta /= two;
        }

        // copy
        for edge_id in 0..graph.num_edges() {
            let i = self.csr.edge_index_to_inside_edge_index[edge_id];
            graph.edges[edge_id].flow = self.csr.inside_edge_list[i].flow;
        }

        Status::Optimal
    }

    fn bfs(&mut self, source: usize, sink: usize, delta: Flow) {
        self.que.clear();
        self.que.push_back(sink);
        let n = self.csr.num_nodes;
        self.csr.distances.fill(n);
        self.csr.distances[sink] = 0;

        while let Some(v) = self.que.pop_front() {
            for e in self.csr.inside_edge_list[self.csr.start[v]..self.csr.start[v + 1]].iter() {
                // e.to -> v
                if self.csr.inside_edge_list[e.rev].residual_capacity() >= delta && self.csr.distances[e.to] == self.csr.num_nodes {
                    self.csr.distances[e.to] = self.csr.distances[v] + 1;
                    if e.to != source {
                        self.que.push_back(e.to);
                    }
                }
            }
        }
    }

    fn dfs(&mut self, u: usize, sink: usize, upper: Flow, delta: Flow) -> Option<Flow> {
        if u == sink {
            return Some(upper);
        }

        let mut res = Flow::zero();
        for i in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = i;
            let v = self.csr.inside_edge_list[i].to;
            let residual_capacity = self.csr.inside_edge_list[i].residual_capacity();

            if !self.csr.is_admissible_edge(u, i) || residual_capacity < delta {
                continue;
            }

            if let Some(d) = self.dfs(v, sink, residual_capacity.min(upper - res), delta) {
                self.csr.push_flow(i, d);
                res += d;
                if res == upper {
                    return Some(res);
                }
            }
        }
        self.current_edge[u] = self.csr.start[u + 1];
        self.csr.distances[u] = self.csr.num_nodes;

        Some(res)
    }
}
