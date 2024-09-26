use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::spanning_tree_structure::{EdgeState, SpanningTreeStructure};
use crate::minimum_cost_flow::status::Status;
use num_traits::NumAssign;
use std::collections::VecDeque;
use std::ops::Neg;

#[derive(Default)]
pub struct ParametricNetworkSimplex<Flow> {
    st: SpanningTreeStructure<Flow>,
    sink: usize,
}

impl<Flow> ParametricNetworkSimplex<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy + Default,
{
    pub fn solve(&mut self, graph: &mut Graph<Flow>) -> Status {
        if graph.is_unbalance() {
            return Status::Unbalanced;
        }

        let (source, sink, artificial_nodes, artificial_edges) = graph.construct_extend_network_one_supply_one_demand();
        self.st.build(graph);
        (self.st.root, self.sink) = (source, sink);

        if !self.make_initial_spanning_tree_structure() {
            // there is no s-t path
            let status = if self.st.satisfy_constraints() { Status::Optimal } else { Status::Infeasible };
            graph.remove_artificial_sub_graph(&artificial_nodes, &artificial_edges);
            return status;
        }
        debug_assert!(self.st.satisfy_optimality_conditions());

        self.run();

        let status = if self.st.satisfy_constraints() { Status::Optimal } else { Status::Infeasible };
        // copy
        graph.excesses = self.st.excesses.clone();
        for edge_id in 0..graph.num_edges() {
            graph.edges[edge_id].flow = self.st.edges[edge_id].flow;
        }
        graph.remove_artificial_sub_graph(&artificial_nodes, &artificial_edges);
        status
    }

    pub(crate) fn run(&mut self) {
        while let Some((leaving_edge_id, delta)) = self.select_leaving_edge() {
            let leaving_edge = &self.st.edges[leaving_edge_id];
            let t2_now_root = if self.st.nodes[leaving_edge.from].parent == leaving_edge.to {
                leaving_edge.from
            } else {
                leaving_edge.to
            };

            self.st.update_flow_in_path(self.st.root, self.sink, delta);
            if self.st.excesses[self.st.root] == Flow::zero() {
                break;
            }

            if let Some((entering_edge_id, t2_new_root)) = self.select_entering_edge_id(leaving_edge_id, t2_now_root) {
                self.dual_pivot(leaving_edge_id, entering_edge_id, t2_now_root, t2_new_root);
                debug_assert!(self.st.satisfy_optimality_conditions());
            } else {
                break;
            }
        }
    }

    // T: shortest path
    // L: A \ T
    // U: empty
    fn make_initial_spanning_tree_structure(&mut self) -> bool {
        let (distances, prev_edge_id) = self.st.shortest_path(self.st.root);

        // there is no s-t path
        if prev_edge_id[self.sink].is_none() {
            return false;
        }

        // make tree structure
        let mut children = vec![Vec::new(); self.st.num_nodes];
        for edge_id in prev_edge_id.iter().filter_map(|&edge_id| edge_id) {
            let edge = &mut self.st.edges[edge_id];
            edge.state = EdgeState::Tree;
            (self.st.nodes[edge.to].parent, self.st.nodes[edge.to].parent_edge_id) = (edge.from, edge_id);
            children[edge.from].push(edge.to);
        }
        (self.st.nodes[self.st.root].parent, self.st.nodes[self.st.root].parent_edge_id) = (usize::MAX, usize::MAX);
        self.st.last_descendent_dft = (0..self.st.num_nodes).collect();

        let mut prev_node = usize::MAX;
        let mut stack = VecDeque::from([(self.st.root, usize::MAX)]);
        let mut seen = vec![false; self.st.num_nodes];
        while let Some((u, parent)) = stack.pop_back() {
            if seen[u] {
                self.st.num_successors[u] += 1;
                if parent != usize::MAX {
                    self.st.last_descendent_dft[parent] = self.st.last_descendent_dft[u];
                    self.st.num_successors[self.st.nodes[u].parent] += self.st.num_successors[u];
                }
                continue;
            }

            seen[u] = true;
            self.st.prev_node_dft[u] = prev_node;
            if prev_node != usize::MAX {
                self.st.next_node_dft[prev_node] = u;
            }
            prev_node = u;
            stack.push_back((u, parent));
            for &child in children[u].iter().rev() {
                stack.push_back((child, u));
            }
        }
        self.st.next_node_dft[prev_node] = self.st.root;

        // determine potentials
        for (u, node) in self.st.nodes.iter_mut().enumerate() {
            node.potential = -distances[u];
        }

        true
    }

    fn select_leaving_edge(&self) -> Option<(usize, Flow)> {
        let mut leaving_edge_id = None;
        let mut mini_delta = Flow::zero();
        let mut now = self.sink;
        while now != self.st.root {
            let (parent, edge_id) = (self.st.nodes[now].parent, self.st.nodes[now].parent_edge_id);
            let edge = &self.st.edges[edge_id];
            assert_eq!(edge.state, EdgeState::Tree);

            let delta = if edge.from == parent { edge.residual_capacity() } else { edge.flow };
            // select the edge closest to the source as the leaving edge
            if leaving_edge_id.is_none() || delta <= mini_delta {
                mini_delta = delta;
                leaving_edge_id = Some(edge_id);
            }

            now = parent;
        }
        Some((leaving_edge_id?, mini_delta.min(self.st.excesses[self.st.root])))
    }

    fn select_entering_edge_id(&self, leaving_edge_id: usize, t2_now_root: usize) -> Option<(usize, usize)> {
        let mut is_t1_node = vec![false; self.st.num_nodes];
        let mut now = self.st.root;
        loop {
            is_t1_node[now] = true;
            now = self.st.next_node_dft[now];
            if now == t2_now_root {
                now = self.st.next_node_dft[self.st.last_descendent_dft[now]];
            }
            if now == self.st.root {
                break;
            }
        }

        let mut entering_edge_id = None;
        let mut t2_new_root = None;
        let mut mini_delta = Flow::zero();
        for (edge_id, edge) in self.st.edges.iter().enumerate() {
            if edge_id == leaving_edge_id {
                continue;
            }

            // t1 -> t2 and lower
            if is_t1_node[edge.from] && !is_t1_node[edge.to] && edge.state == EdgeState::Lower && edge.upper != Flow::zero() {
                let reduced_cost = self.st.reduced_cost(edge);
                if reduced_cost < mini_delta || entering_edge_id.is_none() {
                    mini_delta = reduced_cost;
                    entering_edge_id = Some(edge_id);
                    t2_new_root = Some(edge.to);
                }
            }

            // t2 -> t1 and upper
            if !is_t1_node[edge.from] && is_t1_node[edge.to] && edge.state == EdgeState::Upper {
                let reduced_cost = -self.st.reduced_cost(edge);
                if reduced_cost < mini_delta || entering_edge_id.is_none() {
                    mini_delta = reduced_cost;
                    entering_edge_id = Some(edge_id);
                    t2_new_root = Some(edge.from);
                }
            }
        }
        Some((entering_edge_id?, t2_new_root?))
    }

    fn dual_pivot(&mut self, leaving_edge_id: usize, entering_edge_id: usize, t2_now_root: usize, t2_new_root: usize) {
        if leaving_edge_id == entering_edge_id {
            self.st.edges[entering_edge_id].state = match self.st.edges[entering_edge_id].state {
                EdgeState::Upper => EdgeState::Lower,
                EdgeState::Lower => EdgeState::Upper,
                _ => panic!("state of entering edge {entering_edge_id} is invalid."),
            };
            return;
        }

        // drop leaving edge
        self.st.detach_tree(self.st.root, t2_now_root, leaving_edge_id);

        // enter entering edge
        let entering_edge = &mut self.st.edges[entering_edge_id];
        let attach_node = entering_edge.opposite_side(t2_new_root);
        self.st.re_rooting(t2_now_root, t2_new_root, entering_edge_id);
        self.st.attach_tree(self.st.root, attach_node, t2_new_root, entering_edge_id);
    }
}
