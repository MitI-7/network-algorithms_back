use num_traits::NumAssign;
use std::collections::VecDeque;
use std::fmt::Debug;

#[derive(PartialEq, Debug, Clone)]
pub struct Edge<Flow> {
    pub from: usize,
    pub to: usize,
    pub flow: Flow,
    pub upper: Flow,
}

#[derive(Default)]
pub struct Graph<Flow> {
    num_nodes: usize,
    num_edges: usize,
    pub(crate) edges: Vec<Edge<Flow>>,
    pub(crate) excesses: Vec<Flow>,
}

impl<Flow> Graph<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    #[inline]
    pub fn num_nodes(&self) -> usize {
        self.num_nodes
    }

    #[inline]
    pub fn num_edges(&self) -> usize {
        self.num_edges
    }

    pub fn add_node(&mut self) -> usize {
        self.excesses.push(Flow::zero());
        self.num_nodes += 1;
        self.num_nodes - 1
    }

    pub fn add_nodes(&mut self, num_nodes: usize) -> Vec<usize> {
        self.excesses.extend(vec![Flow::zero(); num_nodes]);
        self.num_nodes += num_nodes;
        ((self.num_nodes - num_nodes)..self.num_nodes).collect()
    }

    // return edge index
    pub fn add_directed_edge(&mut self, from: usize, to: usize, upper: Flow) -> Option<usize> {
        if from >= self.num_nodes || to >= self.num_nodes {
            return None;
        }

        self.edges.push(Edge { from, to, flow: Flow::zero(), upper });

        self.num_edges += 1;
        Some(self.num_edges - 1)
    }

    pub fn get_edge(&self, edge_id: usize) -> Option<Edge<Flow>> {
        if edge_id >= self.edges.len() {
            return None;
        }
        let edge = &self.edges[edge_id];
        Some(Edge { from: edge.from, to: edge.to, flow: edge.flow, upper: edge.upper })
    }

    pub fn maximum_flow(&self, source: usize) -> Flow {
        (0..self.num_edges).fold(Flow::zero(), |mut flow, edge_index| {
            let edge = self.get_edge(edge_index).unwrap();
            if edge.from == source {
                flow += edge.flow;
            } else if edge.to == source {
                flow -= edge.flow;
            }
            flow
        })
    }

    pub fn minimum_cut(&self, source: usize) -> Vec<usize> {
        let mut cut = Vec::new();
        let mut visited = vec![false; self.num_nodes];
        let mut que = VecDeque::from([source]);

        if self.num_edges == 0 {
            return vec![source];
        }

        while let Some(u) = que.pop_front() {
            cut.push(u);
            visited[u] = true;

            // for e in self.neighbors(u) {
            //     if !visited[e.to] && e.residual_capacity() != Flow::zero() {
            //         que.push_back(e.to);
            //     }
            // }
        }

        cut
    }
}
