use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;
use std::collections::VecDeque;

#[derive(Default)]
pub struct EdmondsKarp<Flow> {
    csr: CSR<Flow>,
}

impl<Flow> EdmondsKarp<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        self.csr.build(graph);
        let mut prev = vec![(usize::MAX, usize::MAX); self.csr.num_nodes];
        let mut visited = vec![false; self.csr.num_nodes];

        loop {
            prev.fill((usize::MAX, usize::MAX));
            visited.fill(false);

            // bfs
            let mut queue = VecDeque::from([source]);
            while let Some(u) = queue.pop_front() {
                visited[u] = true;
                if u == sink {
                    break;
                }

                for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
                    let edge = &self.csr.inside_edge_list[edge_id];
                    if visited[edge.to] || edge.residual_capacity() == Flow::zero() {
                        continue;
                    }

                    queue.push_back(edge.to);
                    prev[edge.to] = (u, edge_id);
                }
            }

            if !visited[sink] {
                break;
            }

            // calculate delta
            let mut delta = self.csr.inside_edge_list[prev[sink].1].residual_capacity();
            let mut v = sink;
            while v != source {
                let (u, edge_id) = prev[v];
                delta = delta.min(self.csr.inside_edge_list[edge_id].residual_capacity());
                v = u;
            }

            // update flow
            let mut v = sink;
            while v != source {
                let (u, edge_id) = prev[v];
                self.csr.push_flow(edge_id, delta);
                v = u;
            }
        }

        self.csr.set_flow(graph);
        Status::Optimal
    }
}
