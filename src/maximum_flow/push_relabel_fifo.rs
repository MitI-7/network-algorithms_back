use crate::maximum_flow::csr::CSR;
use crate::maximum_flow::graph::Graph;
use crate::maximum_flow::status::Status;
use num_traits::NumAssign;
use std::collections::VecDeque;

#[derive(Default)]
pub struct PushRelabelFIFO<Flow> {
    csr: CSR<Flow>,
    excesses: Vec<Flow>,

    alpha: usize,
    relabel_count: usize,
    active_nodes: VecDeque<usize>,
    current_edge: Vec<usize>,
    distance_count: Vec<usize>,
}

impl<Flow> PushRelabelFIFO<Flow>
where
    Flow: NumAssign + Ord + Copy + Default,
{
    pub fn new(&mut self, alpha: usize) -> Self {
        Self { csr: CSR::default(), excesses: Vec::new(), alpha, relabel_count: 0, active_nodes: VecDeque::new(), current_edge: Vec::new(), distance_count: Vec::new() }
    }

    pub fn solve(&mut self, source: usize, sink: usize, graph: &mut Graph<Flow>) -> Status {
        if source >= graph.num_nodes() || sink >= graph.num_nodes() || source == sink {
            return Status::BadInput;
        }
        self.csr.build(graph);

        self.pre_process(source, sink);

        while let Some(u) = self.active_nodes.pop_front() {
            // no path to sink
            if u == source || u == sink || self.csr.distances[u] >= self.csr.num_nodes {
                continue;
            }
            self.discharge(u);

            if self.alpha != 0 && self.relabel_count > self.alpha * self.csr.num_nodes {
                self.relabel_count = 0;
                self.csr.update_distances(source, sink);
            }
        }

        self.push_flow_excess_back_to_source(source, sink);

        self.csr.set_flow(graph);

        Status::Optimal
    }

    fn pre_process(&mut self, source: usize, sink: usize) {
        self.excesses.resize(self.csr.num_nodes, Flow::zero());
        self.current_edge.resize(self.csr.num_nodes, 0);
        self.distance_count.resize(self.csr.num_nodes + 1, 0);

        self.csr.update_distances(source, sink);
        self.csr.distances[source] = self.csr.num_nodes;

        for u in 0..self.csr.num_nodes {
            self.distance_count[self.csr.distances[u]] += 1;
            self.current_edge[u] = self.csr.start[u];
        }

        for inside_edge_index in self.csr.start[source]..self.csr.start[source + 1] {
            let edge = &self.csr.inside_edge_list[inside_edge_index];
            let delta = edge.residual_capacity();
            self.excesses[edge.to] += delta;
            self.csr.push_flow(inside_edge_index, delta);
        }

        for u in 0..self.csr.num_nodes {
            if u != source && u != sink && self.excesses[u] > Flow::zero() {
                self.active_nodes.push_back(u);
            }
        }
    }

    fn discharge(&mut self, u: usize) {
        // push
        for edge_id in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = edge_id;
            if self.excesses[u] > Flow::zero() {
                self.push(u, edge_id);
            }

            if self.excesses[u] == Flow::zero() {
                return;
            }
        }
        self.current_edge[u] = self.csr.start[u];

        // relabel
        if self.distance_count[self.csr.distances[u]] == 1 {
            self.gap_relabeling(self.csr.distances[u]);
        } else {
            self.relabel(u);
        }

        if self.excesses[u] > Flow::zero() {
            self.active_nodes.push_back(u);
        }
    }

    // push from u
    fn push(&mut self, u: usize, edge_id: usize) {
        let to = self.csr.inside_edge_list[edge_id].to;
        let delta = self.excesses[u].min(self.csr.inside_edge_list[edge_id].residual_capacity());
        if self.csr.is_admissible_edge(u, edge_id) && delta > Flow::zero() {
            self.csr.push_flow(edge_id, delta);
            self.excesses[u] -= delta;
            self.excesses[to] += delta;
            if self.excesses[to] == delta {
                self.active_nodes.push_back(to);
            }
        }
    }

    fn relabel(&mut self, u: usize) {
        self.relabel_count += 1;
        self.distance_count[self.csr.distances[u]] -= 1;

        let new_distance = self
            .csr
            .neighbors(u)
            .filter(|edge| edge.residual_capacity() > Flow::zero())
            .map(|edge| self.csr.distances[edge.to] + 1)
            .min()
            .unwrap()
            .min(self.csr.num_nodes);

        // assert!(new_distance > self.graph.distances[u]);
        self.csr.distances[u] = new_distance;
        self.distance_count[self.csr.distances[u]] += 1;
    }

    // gap relabeling heuristic
    // set distance[u] >= k to distance[u] = n
    // O(n)
    fn gap_relabeling(&mut self, k: usize) {
        for u in 0..self.csr.num_nodes {
            if self.csr.distances[u] >= k {
                self.distance_count[self.csr.distances[u]] -= 1;
                self.csr.distances[u] = self.csr.distances[u].max(self.csr.num_nodes);
                self.distance_count[self.csr.distances[u]] += 1;
            }
        }
    }

    fn push_flow_excess_back_to_source(&mut self, source: usize, sink: usize) {
        for u in 0..self.csr.num_nodes {
            if u == source || u == sink {
                continue;
            }
            while self.excesses[u] > Flow::zero() {
                let mut visited = vec![false; self.csr.num_nodes];
                self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
                let d = self.dfs(u, source, self.excesses[u], &mut visited);
                self.excesses[u] -= d;
                self.excesses[source] += d;
            }
        }
    }

    fn dfs(&mut self, u: usize, source: usize, flow: Flow, visited: &mut Vec<bool>) -> Flow {
        if u == source {
            return flow;
        }
        visited[u] = true;

        for i in self.current_edge[u]..self.csr.start[u + 1] {
            self.current_edge[u] = i;
            let to = self.csr.inside_edge_list[i].to;
            let residual_capacity = self.csr.inside_edge_list[i].residual_capacity();
            if visited[to] || residual_capacity == Flow::zero() {
                continue;
            }

            let delta = self.dfs(to, source, flow.min(residual_capacity), visited);
            if delta > Flow::zero() {
                self.csr.push_flow(i, delta);
                return delta;
            }
        }
        Flow::zero()
    }
}
