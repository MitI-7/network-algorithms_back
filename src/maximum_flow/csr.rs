use crate::maximum_flow::graph::Graph;
use std::collections::VecDeque;
use std::fmt::Debug;
use std::ops::Sub;

use num_traits::NumAssign;

#[derive(Default, PartialEq, Debug)]
pub struct InsideEdge<Flow> {
    pub to: usize,
    pub flow: Flow,
    pub upper: Flow,
    pub rev: usize,
}

impl<Flow> InsideEdge<Flow>
where
    Flow: Sub<Output = Flow> + Copy,
{
    pub fn residual_capacity(&self) -> Flow {
        self.upper - self.flow
    }
}

#[derive(Default)]
pub struct CSR<Flow> {
    pub num_nodes: usize,
    pub num_edges: usize,
    pub edge_index_to_inside_edge_index: Vec<usize>,

    pub start: Vec<usize>,
    pub inside_edge_list: Vec<InsideEdge<Flow>>,
    pub distances: Vec<usize>, // distance from u to sink in residual network
    que: VecDeque<usize>,
}

impl<Flow> CSR<Flow>
where
    Flow: NumAssign + Ord + Copy,
{
    pub fn build(&mut self, graph: &mut Graph<Flow>) {
        self.num_nodes = graph.num_nodes();
        self.num_edges = graph.num_edges();

        // initialize
        self.edge_index_to_inside_edge_index.resize(self.num_edges, usize::MAX);
        self.start.resize(self.num_nodes + 1, 0);
        self.inside_edge_list = (0..2 * self.num_edges).map(|_| InsideEdge { to: 0, flow: Flow::zero(), upper: Flow::zero(), rev: 0 }).collect();
        self.distances.resize(self.num_nodes, self.num_nodes);

        let mut degree = vec![0; self.num_nodes];
        for edge in graph.edges.iter() {
            degree[edge.to] += 1;
            degree[edge.from] += 1;
        }

        for i in 1..=self.num_nodes {
            self.start[i] += self.start[i - 1] + degree[i - 1];
        }

        let mut counter = vec![0; self.num_nodes];
        for (edge_index, e) in graph.edges.iter().enumerate() {
            let (u, v) = (e.from, e.to);
            let inside_edge_index_u = self.start[u] + counter[u];
            counter[u] += 1;
            let inside_edge_index_v = self.start[v] + counter[v];
            self.edge_index_to_inside_edge_index[edge_index] = inside_edge_index_u;
            counter[v] += 1;

            self.inside_edge_list[inside_edge_index_u] = InsideEdge { to: v, flow: Flow::zero(), upper: e.upper, rev: inside_edge_index_v };
            self.inside_edge_list[inside_edge_index_v] = InsideEdge { to: u, flow: e.upper, upper: e.upper, rev: inside_edge_index_u };
        }
    }

    pub fn set_flow(&self, graph: &mut Graph<Flow>) {
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
    pub fn push_flow(&mut self, inside_edge_index: usize, flow: Flow) {
        let rev = self.inside_edge_list[inside_edge_index].rev;

        // update flow
        self.inside_edge_list[inside_edge_index].flow += flow;
        self.inside_edge_list[rev].flow -= flow;
    }

    // O(n + m)
    // calculate the distance from u to sink in the residual network
    // if such a path does not exist, distance[u] becomes self.num_nodes
    pub fn update_distances(&mut self, source: usize, sink: usize) {
        self.que.clear();
        self.que.push_back(sink);
        self.distances.fill(self.num_nodes);
        self.distances[sink] = 0;

        while let Some(v) = self.que.pop_front() {
            for e in self.inside_edge_list[self.start[v]..self.start[v + 1]].iter() {
                // e.to -> v
                if e.flow > Flow::zero() && self.distances[e.to] == self.num_nodes {
                    self.distances[e.to] = self.distances[v] + 1;
                    if e.to != source {
                        self.que.push_back(e.to);
                    }
                }
            }
        }
    }

    #[inline]
    pub fn is_admissible_edge(&self, from: usize, i: usize) -> bool {
        self.inside_edge_list[i].residual_capacity() > Flow::zero() && self.distances[from] == self.distances[self.inside_edge_list[i].to] + 1
    }
}
