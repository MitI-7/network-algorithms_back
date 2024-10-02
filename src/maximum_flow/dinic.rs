use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;

#[derive(Default)]
pub struct Dinic<Flow> {
    pub csr: CSR<Flow>,
    current_edge: Vec<usize>,
}

impl<Flow> Dinic<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        self.csr.build(graph);
        self.current_edge.resize(graph.num_nodes(), 0);

        let upper = self.csr.neighbors(source).fold(Flow::zero(), |sum, e| sum + e.upper);
        let mut flow = Flow::zero();
        while flow < upper {
            self.csr.update_distances(source, sink);

            // no s-t path
            if self.csr.distances[source] >= self.csr.num_nodes {
                break;
            }

            self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
            match self.dfs(source, sink, upper) {
                Some(delta) => flow += delta,
                None => break,
            }
        }

        self.csr.set_flow(graph);
        Status::Optimal
    }

    fn dfs(&mut self, u: usize, sink: usize, upper: Flow) -> Option<Flow> {
        if u == sink {
            return Some(upper);
        }

        let mut res = Flow::zero();
        for i in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = i;
            let v = self.csr.inside_edge_list[i].to;
            let residual_capacity = self.csr.inside_edge_list[i].residual_capacity();

            if !self.csr.is_admissible_edge(u, i) {
                continue;
            }

            if let Some(d) = self.dfs(v, sink, residual_capacity.min(upper - res)) {
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
