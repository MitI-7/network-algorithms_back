use num_traits::Float;
use num_traits::ToPrimitive;
use std::fmt::Debug;

#[derive(PartialEq, Debug, Clone)]
pub struct Edge<Flow> {
    pub from: usize,
    pub to: usize,
    pub flow: Flow,
    pub upper: Flow,
    pub gain: Flow,
}

#[derive(Default)]
pub struct Graph<Flow> {
    num_nodes: usize,
    num_edges: usize,
    pub(crate) edges: Vec<Edge<Flow>>,
    b: Vec<Flow>,
    pub(crate) excesses: Vec<Flow>,
}

impl<Flow> Graph<Flow>
where
    Flow: Float + PartialOrd + Copy + Clone + ToPrimitive,
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
        self.b.push(Flow::zero());
        self.excesses.push(Flow::zero());
        self.num_nodes += 1;
        self.num_nodes - 1
    }

    pub fn add_nodes(&mut self, num_nodes: usize) -> Vec<usize> {
        self.b.extend(vec![Flow::zero(); num_nodes]);
        self.excesses.extend(vec![Flow::zero(); num_nodes]);
        self.num_nodes += num_nodes;
        ((self.num_nodes - num_nodes)..self.num_nodes).collect()
    }

    pub fn add_supply(&mut self, u: usize, supply: Flow) {
        self.b[u] = self.b[u] + supply;
        self.excesses[u] = self.excesses[u] + supply;
    }

    pub fn add_demand(&mut self, u: usize, demand: Flow) {
        self.b[u] = self.b[u] - demand;
        self.excesses[u] = self.excesses[u] - demand;
    }

    // return edge index
    pub fn add_directed_edge(&mut self, from: usize, to: usize, upper: Flow, gain: Flow) -> Option<usize> {
        if upper <= Flow::zero() || from >= self.num_nodes || to >= self.num_nodes || gain <= Flow::zero() {
            return None;
        }

        self.edges.push(Edge { from, to, flow: Flow::zero(), upper, gain });

        self.num_edges += 1;
        Some(self.num_edges - 1)
    }

    pub fn get_edge(&self, edge_id: usize) -> Option<Edge<Flow>> {
        if edge_id >= self.edges.len() {
            return None;
        }
        let edge = &self.edges[edge_id];
        Some(Edge { from: edge.from, to: edge.to, flow: edge.flow, upper: edge.upper, gain: edge.gain })
    }

    pub fn maximum_flow(&self, sink: usize) -> Flow {
        (0..self.num_edges()).fold(Flow::zero(), |flow, edge_index| {
            let e = &self.get_edge(edge_index).unwrap();
            flow + if e.to == sink { e.flow * e.gain } else { Flow::zero() }
        })
    }
}
