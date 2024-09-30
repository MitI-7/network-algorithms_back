use num_traits::NumAssign;
use std::fmt::Debug;
use std::ops::Neg;

#[derive(PartialEq, Debug, Clone)]
pub struct Edge<Flow> {
    pub from: usize,
    pub to: usize,
    pub flow: Flow,
    pub lower: Flow,
    pub upper: Flow,
    pub cost: Flow,
}

#[derive(Default)]
pub struct Graph<Flow> {
    num_nodes: usize,
    num_edges: usize,
    pub(crate) edges: Vec<Edge<Flow>>,
    pub(crate) b: Vec<Flow>,
    pub(crate) lowers: Vec<Flow>,
    pub(crate) excesses: Vec<Flow>,
    pub(crate) is_reversed: Vec<bool>,
}

impl<Flow> Graph<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
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
        self.b[u] += supply;
        self.excesses[u] += supply;
    }

    pub fn add_demand(&mut self, u: usize, demand: Flow) {
        self.b[u] -= demand;
        self.excesses[u] -= demand;
    }

    // return edge index
    pub fn add_directed_edge(&mut self, from: usize, to: usize, lower: Flow, upper: Flow, cost: Flow) -> Option<usize> {
        if lower > upper || from >= self.num_nodes || to >= self.num_nodes {
            return None;
        }

        if cost >= Flow::zero() {
            self.edges.push(Edge { from, to, flow: Flow::zero(), lower: Flow::zero(), upper: upper - lower, cost });
            self.excesses[from] -= lower;
            self.excesses[to] += lower;
            self.lowers.push(lower);
            self.is_reversed.push(false);
        } else {
            self.edges.push(Edge { from: to, to: from, flow: Flow::zero(), lower: Flow::zero(), upper: upper - lower, cost: -cost });
            self.excesses[from] -= upper;
            self.excesses[to] += upper;
            self.lowers.push(lower);
            self.is_reversed.push(true);
        }

        self.num_edges += 1;
        Some(self.num_edges - 1)
    }

    pub fn get_edge(&self, edge_id: usize) -> Option<Edge<Flow>> {
        if edge_id >= self.edges.len() {
            return None;
        }
        let edge = &self.edges[edge_id];
        let lower = self.lowers[edge_id];
        if self.is_reversed[edge_id] {
            Some(Edge { from: edge.to, to: edge.from, flow: edge.upper - edge.flow + lower, lower, upper: edge.upper + lower, cost: -edge.cost })
        } else {
            Some(Edge { from: edge.from, to: edge.to, flow: edge.flow + lower, lower, upper: edge.upper + lower, cost: edge.cost })
        }
    }

    pub fn minimum_cost(&self) -> Flow {
        (0..self.num_edges).fold(Flow::zero(), |cost, edge_id| {
            let edge = self.get_edge(edge_id).unwrap();
            cost + edge.cost * edge.flow
        })
    }

    pub fn is_unbalance(&self) -> bool {
        self.b.iter().fold(Flow::zero(), |sum, &excess| sum + excess) != Flow::zero()
    }

    pub(crate) fn construct_extend_network_one_supply_one_demand(&mut self) -> (usize, usize, Vec<usize>, Vec<usize>) {
        let mut artificial_edges = Vec::new();
        let (source, sink) = (self.add_node(), self.add_node());
        for u in 0..self.num_nodes() {
            if u == source || u == sink {
                continue;
            }
            if self.excesses[u] > Flow::zero() {
                artificial_edges.push(self.add_directed_edge(source, u, Flow::zero(), self.excesses[u], Flow::zero()).unwrap());
                self.excesses[source] = self.excesses[source] + self.excesses[u];
            }
            if self.excesses[u] < Flow::zero() {
                artificial_edges.push(self.add_directed_edge(u, sink, Flow::zero(), -self.excesses[u], Flow::zero()).unwrap());
                self.excesses[sink] = self.excesses[sink] + self.excesses[u];
            }
            self.excesses[u] = Flow::zero();
        }

        (source, sink, vec![source, sink], artificial_edges)
    }

    pub(crate) fn construct_extend_network_feasible_solution(&mut self) -> (usize, Vec<usize>, Vec<usize>) {
        let inf_cost = self.edges.iter().map(|e| e.cost).fold(Flow::one(), |acc, cost| acc + cost); // all edge costs are non-negative

        // add artificial nodes
        let root = self.add_node();

        // add artificial edges
        let mut artificial_edges = Vec::new();
        for u in 0..self.num_nodes {
            if u == root {
                continue;
            }

            let excess = self.excesses[u];
            if excess >= Flow::zero() {
                // u -> root
                let edge_id = self.add_directed_edge(u, root, Flow::zero(), excess, inf_cost).unwrap();
                self.edges[edge_id].flow = excess;
                artificial_edges.push(edge_id);
            } else {
                // root -> u
                let edge_id = self.add_directed_edge(root, u, Flow::zero(), -excess, inf_cost).unwrap();
                self.edges[edge_id].flow = -excess;
                artificial_edges.push(edge_id);
            }
            self.excesses[u] = Flow::zero();
        }

        (root, vec![root], artificial_edges)
    }

    pub(crate) fn remove_artificial_sub_graph(&mut self, artificial_nodes: &[usize], artificial_edges: &[usize]) {
        self.edges.truncate(self.num_edges - artificial_edges.len());
        self.b.truncate(self.num_nodes - artificial_nodes.len());
        self.lowers.truncate(self.num_edges - artificial_edges.len());
        self.excesses.truncate(self.num_nodes - artificial_nodes.len());
        self.is_reversed.truncate(self.num_edges - artificial_edges.len());

        self.num_nodes -= artificial_nodes.len();
        self.num_edges -= artificial_edges.len();
    }
}
