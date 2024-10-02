use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;

#[derive(Default)]
pub struct FordFulkerson<Flow> {
    csr: CSR<Flow>,
}

impl<Flow> FordFulkerson<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        self.csr.build(graph);
        let mut visited = vec![false; self.csr.num_nodes];

        let upper = self.csr.neighbors(source).fold(Flow::zero(), |sum, e| sum + e.upper);
        let mut flow = Flow::zero();
        loop {
            visited.fill(false);
            match self.dfs(source, sink, upper, &mut visited) {
                Some(delta) => flow += delta,
                None => break,
            }
        }

        self.csr.set_flow(graph);
        Status::Optimal
    }

    fn dfs(&mut self, u: usize, sink: usize, flow: Flow, visited: &mut Vec<bool>) -> Option<Flow> {
        if u == sink {
            return Some(flow);
        }
        visited[u] = true;

        for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
            let edge = &self.csr.inside_edge_list[edge_id];
            if visited[edge.to] || edge.residual_capacity() == Flow::zero() {
                continue;
            }

            if let Some(d) = self.dfs(edge.to, sink, flow.min(edge.residual_capacity()), visited) {
                self.csr.push_flow(edge_id, d);
                return Some(d);
            }
        }
        None
    }
}
