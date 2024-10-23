use crate::maximum_flow::capacity_scaling::CapacityScaling;
use crate::maximum_flow::graph;
use crate::minimum_cost_flow::csr::{InsideEdge, CSR};
use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::status::Status;
use num_traits::{FromPrimitive, NumAssign};
use std::collections::VecDeque;
use std::ops::Neg;

pub struct CostScalingPushRelabel<Flow> {
    csr: CSR<Flow>,
    active_nodes: VecDeque<usize>,
    current_edge: Vec<usize>,
    alpha: Flow,
}

impl<Flow> Default for CostScalingPushRelabel<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy + FromPrimitive + Default,
{
    fn default() -> Self {
        Self { csr: CSR::default(), active_nodes: VecDeque::new(), current_edge: Vec::new(), alpha: Flow::from_isize(16).unwrap() }
    }
}

#[allow(dead_code)]
impl<Flow> CostScalingPushRelabel<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy + FromPrimitive + Default,
{
    // scaling_factor: it was usually between 8 and 24. default scaling factor is 16
    pub fn new(scaling_factor: Flow) -> Self {
        assert!(scaling_factor > Flow::one());
        Self { csr: CSR::default(), active_nodes: VecDeque::new(), current_edge: Vec::new(), alpha: scaling_factor }
    }

    pub fn solve(&mut self, graph: &mut Graph<Flow>) -> Status {
        if graph.is_unbalance() {
            return Status::Unbalanced;
        }
        self.csr.build(graph);

        // all edge costs are non-negative
        if self.csr.excesses.iter().all(|&excess| excess == Flow::zero()) {
            return Status::Optimal;
        }

        if !self.check_feasibility(graph) {
            return Status::Infeasible;
        }

        self.current_edge.resize(self.csr.num_nodes, 0);
        let gamma = self.csr.inside_edge_list.iter().map(|e| e.cost).max().unwrap_or(Flow::one()); // all edge costs are non-negative
        let cost_scaling_factor = self.alpha * Flow::from_usize(self.csr.num_nodes).unwrap();
        let mut epsilon = Flow::one().max(gamma * cost_scaling_factor);

        // scale cost
        self.csr.inside_edge_list.iter_mut().for_each(|e| e.cost *= cost_scaling_factor);
        loop {
            epsilon = Flow::one().max(epsilon / self.alpha);
            self.refine(epsilon);
            if epsilon == Flow::one() {
                break;
            }
        }
        // unscale cost
        self.csr.inside_edge_list.iter_mut().for_each(|e| e.cost /= cost_scaling_factor);

        self.csr.set_flow(graph);

        Status::Optimal
    }

    // make epsilon-optimal flow
    fn refine(&mut self, epsilon: Flow) {
        // make 0-optimal pseudo flow
        for u in 0..self.csr.num_nodes {
            for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
                let edge = &self.csr.inside_edge_list[edge_id];

                let reduced_cost = self.csr.reduced_cost(u, edge);
                if reduced_cost < Flow::zero() {
                    self.csr.push_flow(u, edge_id, edge.residual_capacity());
                    debug_assert!(self.csr.inside_edge_list[edge_id].flow == self.csr.inside_edge_list[edge_id].upper);
                } else if reduced_cost > Flow::zero() {
                    self.csr.push_flow(u, edge_id, -edge.flow);
                    debug_assert!(self.csr.inside_edge_list[edge_id].flow == Flow::zero());
                }
            }
        }

        self.current_edge.iter_mut().enumerate().for_each(|(u, e)| *e = self.csr.start[u]);
        debug_assert_eq!(self.active_nodes.len(), 0);
        self.active_nodes.extend((0..self.csr.num_nodes).filter(|&u| self.csr.excesses[u] > Flow::zero()));

        // 0-optimal pseudo flow -> epsilon-optimal feasible flow
        while let Some(u) = self.active_nodes.pop_back() {
            self.discharge(u, epsilon);
        }
    }

    fn discharge(&mut self, u: usize, epsilon: Flow) {
        while self.csr.excesses[u] > Flow::zero() {
            self.push(u, epsilon);

            if self.csr.excesses[u] == Flow::zero() {
                break;
            }

            self.relabel(u, epsilon);
        }
    }

    fn is_admissible(&self, u: usize, edge: &InsideEdge<Flow>, _epsilon: Flow) -> bool {
        self.csr.reduced_cost(u, edge) < Flow::zero()
    }

    fn push(&mut self, u: usize, epsilon: Flow) {
        debug_assert!(self.csr.excesses[u] > Flow::zero());

        for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
            let edge = &self.csr.inside_edge_list[edge_id];
            let to = edge.to;
            if edge.residual_capacity() <= Flow::zero() {
                continue;
            }

            if !self.is_admissible(u, edge, epsilon) {
                continue;
            }

            if !self.look_ahead(to, epsilon) {
                if !self.is_admissible(u, &self.csr.inside_edge_list[edge_id], epsilon) {
                    continue;
                }
            }

            let flow = self.csr.inside_edge_list[edge_id].residual_capacity().min(self.csr.excesses[u]);
            self.csr.push_flow(u, edge_id, flow);

            if self.csr.excesses[to] > Flow::zero() && self.csr.excesses[to] <= flow {
                self.active_nodes.push_back(to);
            }

            if self.csr.excesses[u] == Flow::zero() {
                self.current_edge[u] = edge_id;
                return;
            }
        }

        // node has no admissible edge
        self.current_edge[u] = self.csr.start[u];
    }

    fn relabel(&mut self, u: usize, epsilon: Flow) {
        let guaranteed_new_potential = self.csr.potentials[u] + epsilon;

        let mut mini_potential = None;
        let mut previous_mini_potential = None;
        let mut current_edges_for_u = 0;

        for edge_id in self.csr.start[u]..self.csr.start[u + 1] {
            if self.csr.inside_edge_list[edge_id].residual_capacity() <= Flow::zero() {
                continue;
            }

            let to = self.csr.inside_edge_list[edge_id].to;
            let cost = self.csr.inside_edge_list[edge_id].cost;

            let new_potential = self.csr.potentials[to] + cost;
            if mini_potential.is_none() || new_potential < mini_potential.unwrap() {
                // adding epsilon creates an admissible edge
                if new_potential < guaranteed_new_potential {
                    self.csr.potentials[u] = guaranteed_new_potential;
                    self.current_edge[u] = edge_id;
                    return;
                }

                previous_mini_potential = mini_potential;
                mini_potential = Some(new_potential);
                current_edges_for_u = edge_id;
            }
        }

        // increasing the potential cannot create an admissible edge.
        if mini_potential.is_none() {
            if self.csr.excesses[u] != Flow::zero() {
                return;
            } else {
                // the potential can be lowered as much as desired, but guaranteed_new_potential is applied
                self.csr.potentials[u] = guaranteed_new_potential;
                self.current_edge[u] = 0;
            }
            return;
        }

        // adding more than epsilon creates an admissible edge
        let new_potential = mini_potential.unwrap() + epsilon;
        self.csr.potentials[u] = new_potential;

        self.current_edge[u] = if previous_mini_potential.is_none() || previous_mini_potential.unwrap() >= new_potential {
            current_edges_for_u
        } else {
            self.csr.start[u]
        };
    }

    fn look_ahead(&mut self, u: usize, epsilon: Flow) -> bool {
        if self.csr.excesses[u] < Flow::zero() {
            return true;
        }

        // search admissible edge
        for edge_id in self.current_edge[u]..self.csr.start[u + 1] {
            let edge = &self.csr.inside_edge_list[edge_id];
            if edge.residual_capacity() <= Flow::zero() {
                continue;
            }

            if self.is_admissible(u, edge, epsilon) {
                self.current_edge[u] = edge_id;
                return true;
            }
        }

        self.relabel(u, epsilon);
        false
    }

    fn check_feasibility(&self, graph: &Graph<Flow>) -> bool {
        let mut maximum_flow_graph = graph::Graph::default();
        maximum_flow_graph.add_nodes(graph.num_nodes());
        let source = maximum_flow_graph.add_node();
        let sink = maximum_flow_graph.add_node();

        let mut excesses = graph.b.clone();
        for (edge_id, edge) in graph.edges.iter().enumerate() {
            let (from, to) = if graph.is_reversed[edge_id] { (edge.to, edge.from) } else { (edge.from, edge.to) };
            excesses[from] -= graph.lowers[edge_id];
            excesses[to] += graph.lowers[edge_id];
            maximum_flow_graph.add_directed_edge(from, to, edge.upper);
        }

        let mut total_excess = Flow::zero();
        for u in 0..graph.num_nodes() {
            if excesses[u] > Flow::zero() {
                maximum_flow_graph.add_directed_edge(source, u, excesses[u]);
                total_excess += excesses[u];
            }
            if excesses[u] < Flow::zero() {
                maximum_flow_graph.add_directed_edge(u, sink, -excesses[u]);
            }
        }
        CapacityScaling::default().solve(source, sink, &mut maximum_flow_graph);
        maximum_flow_graph.maximum_flow(source) >= total_excess
    }
}
