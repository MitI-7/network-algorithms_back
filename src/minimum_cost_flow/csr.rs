use crate::minimum_cost_flow::graph::Graph;
use num_traits::NumAssign;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fmt::Debug;
use std::ops::Neg;

#[allow(dead_code)]
impl<Flow> InsideEdge<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn residual_capacity(&self) -> Flow {
        self.upper - self.flow
    }

    pub fn is_feasible(&self) -> bool {
        Flow::zero() <= self.flow && self.flow <= self.upper
    }
}

#[derive(Default)]
pub struct CSR<Flow> {
    pub num_nodes: usize,
    pub num_edges: usize,
    pub edge_index_to_inside_edge_index: Vec<usize>,

    pub excesses: Vec<Flow>,
    pub potentials: Vec<Flow>,

    pub start: Vec<usize>,
    pub inside_edge_list: Vec<InsideEdge<Flow>>,
}

#[derive(Default, Debug)]
pub struct InsideEdge<Flow> {
    pub to: usize,
    pub flow: Flow,
    pub upper: Flow,
    pub cost: Flow,
    pub rev: usize,
}

#[allow(dead_code)]
impl<Flow> CSR<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn build(&mut self, graph: &Graph<Flow>) {
        if graph.num_nodes() == 0 {
            return;
        }

        self.num_nodes = graph.num_nodes();
        self.num_edges = graph.num_edges();
        self.excesses = graph.excesses.clone();

        // initialize
        self.edge_index_to_inside_edge_index.resize(self.num_edges, usize::MAX);
        self.start.resize(self.num_nodes + 1, 0);
        self.inside_edge_list = (0..2 * self.num_edges)
            .map(|_| InsideEdge { to: 0, flow: Flow::zero(), upper: Flow::zero(), cost: Flow::zero(), rev: 0 })
            .collect();
        self.potentials.resize(self.num_nodes, Flow::zero());

        let mut degree = vec![0; self.num_nodes];
        for edge in graph.edges.iter() {
            degree[edge.to] += 1;
            degree[edge.from] += 1;
        }

        for i in 1..=self.num_nodes {
            self.start[i] += self.start[i - 1] + degree[i - 1];
        }

        let mut counter = vec![0; self.num_nodes];
        for (edge_index, edge) in graph.edges.iter().enumerate() {
            let (u, v) = (edge.from, edge.to);
            let inside_edge_index_u = self.start[u] + counter[u];
            counter[u] += 1;
            let inside_edge_index_v = self.start[v] + counter[v];
            self.edge_index_to_inside_edge_index[edge_index] = inside_edge_index_u;
            counter[v] += 1;

            assert_ne!(inside_edge_index_u, inside_edge_index_v);

            // u -> v
            self.inside_edge_list[inside_edge_index_u] = InsideEdge { to: v, flow: edge.flow, upper: edge.upper, cost: edge.cost, rev: inside_edge_index_v };
            // v -> u
            self.inside_edge_list[inside_edge_index_v] = InsideEdge { to: u, flow: edge.upper - edge.flow, upper: edge.upper, cost: -edge.cost, rev: inside_edge_index_u };

            assert!(edge.cost >= Flow::zero());
            assert!(edge.upper >= Flow::zero());
        }
    }

    pub fn set_flow(&self, graph: &mut Graph<Flow>) {
        graph.excesses = self.excesses.clone();
        for edge_id in 0..graph.num_edges() {
            let i = self.edge_index_to_inside_edge_index[edge_id];
            graph.edges[edge_id].flow = self.inside_edge_list[i].flow;
        }
    }

    #[inline]
    pub fn neighbors(&self, u: usize) -> std::slice::Iter<InsideEdge<Flow>> {
        self.inside_edge_list[self.start[u]..self.start[u + 1]].iter()
    }

    #[inline]
    pub fn push_flow(&mut self, u: usize, edge_id: usize, flow: Flow) {
        let rev = self.inside_edge_list[edge_id].rev;
        let to = self.inside_edge_list[edge_id].to;
        self.inside_edge_list[edge_id].flow += flow;
        self.inside_edge_list[rev].flow -= flow;
        self.excesses[u] -= flow;
        self.excesses[to] += flow;
    }

    pub fn calculate_distance_from_source(&mut self, source: usize) -> (Vec<Option<Flow>>, Vec<Option<usize>>) {
        let mut prev = vec![None; self.num_nodes];
        let mut bh = BinaryHeap::new();
        let mut dist: Vec<Option<Flow>> = vec![None; self.num_nodes];
        let mut visited = vec![false; self.num_nodes];

        bh.push((Reverse(Flow::zero()), source));
        dist[source] = Some(Flow::zero());

        while let Some((d, u)) = bh.pop() {
            if visited[u] {
                continue;
            }
            visited[u] = true;

            for edge_id in self.start[u]..self.start[u + 1] {
                let edge = &self.inside_edge_list[edge_id];
                if edge.residual_capacity() == Flow::zero() {
                    continue;
                }

                let new_dist = d.0 + self.reduced_cost(u, edge);
                if dist[edge.to].is_none() || dist[edge.to].unwrap() > new_dist {
                    dist[edge.to] = Some(new_dist);
                    prev[edge.to] = Some(edge_id);
                    bh.push((Reverse(new_dist), edge.to));
                }
            }
        }

        (dist, prev)
    }

    #[inline]
    pub fn reduced_cost(&self, u: usize, e: &InsideEdge<Flow>) -> Flow {
        e.cost - self.potentials[u] + self.potentials[e.to]
    }

    #[inline]
    pub fn reduced_cost_rev(&self, u: usize, e: &InsideEdge<Flow>) -> Flow {
        -(e.cost - self.potentials[u] + self.potentials[e.to])
    }
}
