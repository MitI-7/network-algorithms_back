use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;

#[derive(Default)]
pub struct ShortestAugmentingPath<Flow> {
    csr: CSR<Flow>,
    pub current_edge: Vec<usize>,
}

impl<Flow> ShortestAugmentingPath<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        self.csr.build(graph);
        self.csr.update_distances(source, sink);
        self.current_edge.resize(self.csr.num_nodes, 0);

        let mut flow = Flow::zero();
        let upper = self.csr.neighbors(source).fold(Flow::zero(), |sum, e| sum + e.upper);
        while self.csr.distances[source] < self.csr.num_nodes {
            self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
            if let Some(delta) = self.dfs(source, sink, upper) {
                flow += delta;
            }
        }

        self.csr.set_flow(graph);
        Status::Optimal
    }

    fn dfs(&mut self, u: usize, sink: usize, upper: Flow) -> Option<Flow> {
        if u == sink {
            return Some(upper);
        }

        for i in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = i;
            let e = &self.csr.inside_edge_list[i];
            if self.csr.is_admissible_edge(u, i) {
                // advance
                if let Some(delta) = self.dfs(e.to, sink, upper.min(e.residual_capacity())) {
                    self.csr.push_flow(i, delta);
                    return Some(delta);
                }
            }
        }

        // retreat
        self.csr.distances[u] = self.csr.num_nodes;
        for e in self.csr.inside_edge_list[self.csr.start[u]..self.csr.start[u + 1]].iter() {
            if e.residual_capacity() > Flow::zero() {
                self.csr.distances[u] = self.csr.distances[u].min(self.csr.distances[e.to] + 1);
            }
        }

        None
    }
}
