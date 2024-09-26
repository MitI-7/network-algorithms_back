use crate::minimum_cost_flow::graph::Graph;
use crate::minimum_cost_flow::network_simplex_pivot_rules::PivotRule;
use crate::minimum_cost_flow::spanning_tree_structure::{EdgeState, InternalEdge, SpanningTreeStructure};
use crate::minimum_cost_flow::status::Status;
use num_traits::NumAssign;
use std::collections::VecDeque;
use std::ops::Neg;

#[derive(Default)]
pub struct DualNetworkSimplex<Flow> {
    st: SpanningTreeStructure<Flow>,
    sink: usize,
}

impl<Flow> DualNetworkSimplex<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn solve<Pivot: PivotRule<Flow>>(&mut self, pivot: &mut Pivot, graph: &mut Graph<Flow>) -> Status {
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

        self.run(pivot);

        let status = if self.st.satisfy_constraints() { Status::Optimal } else { Status::Infeasible };

        // copy
        graph.excesses = self.st.excesses.clone();
        for edge_id in 0..graph.num_edges() {
            graph.edges[edge_id].flow = self.st.edges[edge_id].flow;
        }
        graph.remove_artificial_sub_graph(&artificial_nodes, &artificial_edges);
        status
    }

    fn run<Pivot: PivotRule<Flow>>(&mut self, pivot: &mut Pivot) {
        while let Some(leaving_edge_id) = pivot.find_entering_edge(&self.st, Self::calculate_violation) {
            let leaving_edge = &self.st.edges[leaving_edge_id];
            let t2_now_root = if self.st.nodes[leaving_edge.from].parent == leaving_edge.to {
                leaving_edge.from
            } else {
                leaving_edge.to
            };

            if let Some((entering_edge_id, t2_new_root)) = self.select_entering_edge_id(leaving_edge_id, t2_now_root) {
                let delta = Self::calculate_violation(&self.st.edges[leaving_edge_id], &self.st);
                let apex = self.find_apex(entering_edge_id);

                // update flow
                self.st.update_flow_in_cycle(entering_edge_id, delta, apex);
                assert!(self.st.edges[leaving_edge_id].is_lower() || self.st.edges[leaving_edge_id].is_upper());

                self.dual_pivot(leaving_edge_id, entering_edge_id, t2_now_root, t2_new_root);
                debug_assert!(self.st.validate_num_successors(self.st.root));
                debug_assert!(self.st.satisfy_optimality_conditions());
            } else {
                break;
            }
        }
    }

    fn calculate_violation(edge: &InternalEdge<Flow>, _: &SpanningTreeStructure<Flow>) -> Flow {
        if edge.flow < Flow::zero() {
            -edge.flow
        } else if edge.flow > edge.upper {
            edge.flow - edge.upper
        } else {
            Flow::zero()
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

        // send flow from source to sink
        self.st.update_flow_in_path(self.st.root, self.sink, self.st.excesses[self.st.root]);
        assert!(self.st.excesses[self.st.root] == Flow::zero() && self.st.excesses[self.sink] == Flow::zero());

        true
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

        let flow_direction_t1_t2 = |edge: &InternalEdge<Flow>| {
            // (t1 -> t2 and lower) or (t2 -> t1 and upper)
            (is_t1_node[edge.from] && !is_t1_node[edge.to] && edge.flow <= Flow::zero()) || (!is_t1_node[edge.from] && is_t1_node[edge.to] && edge.flow >= edge.upper)
        };

        let leaving_edge_flow_direction = flow_direction_t1_t2(&self.st.edges[leaving_edge_id]);

        let mut entering_edge_id = None;
        let mut t2_new_root = None;
        let mut mini_delta = Flow::zero();
        for (edge_id, edge) in self.st.edges.iter().enumerate() {
            if edge.state == EdgeState::Tree || edge.upper == Flow::zero() {
                continue;
            }

            let entering_edge_flow_direction = flow_direction_t1_t2(edge);
            if leaving_edge_flow_direction == entering_edge_flow_direction || is_t1_node[edge.from] == is_t1_node[edge.to] {
                continue;
            }

            let reduced_cost = if edge.state == EdgeState::Lower {
                self.st.reduced_cost(edge)
            } else {
                -self.st.reduced_cost(edge)
            };
            assert!(reduced_cost >= Flow::zero());

            if reduced_cost < mini_delta || entering_edge_id.is_none() {
                mini_delta = reduced_cost;
                entering_edge_id = Some(edge_id);
                t2_new_root = if (entering_edge_flow_direction && edge.state == EdgeState::Lower) || (!entering_edge_flow_direction && edge.state == EdgeState::Upper) {
                    Some(edge.to)
                } else {
                    Some(edge.from)
                };
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

        // drop leaving edge and detach tree
        self.st.detach_tree(self.st.root, t2_now_root, leaving_edge_id);

        // if the size of subtree t2 is larger than that of subtree t1, swap t1 and t2.
        let (t1_new_root, t2_new_root, t2_now_root, new_attach_node) = if self.st.num_successors[t2_now_root] * 2 >= self.st.num_nodes {
            (t2_now_root, self.st.edges[entering_edge_id].opposite_side(t2_new_root), self.st.root, t2_new_root)
        } else {
            (self.st.root, t2_new_root, t2_now_root, self.st.edges[entering_edge_id].opposite_side(t2_new_root))
        };

        // enter entering edge and attach tree
        self.st.re_rooting(t2_now_root, t2_new_root, entering_edge_id);
        self.st.attach_tree(t1_new_root, new_attach_node, t2_new_root, entering_edge_id);
        self.st.root = t1_new_root;
        self.st.nodes[self.st.root].parent = usize::MAX;
        assert_eq!(self.st.nodes[self.st.root].parent_edge_id, usize::MAX);
    }

    fn find_apex(&self, entering_edge_id: usize) -> usize {
        let entering_edge = &self.st.edges[entering_edge_id];
        let (mut u, mut v) = (entering_edge.from, entering_edge.to);
        while u != v {
            let (u_num, v_num) = (self.st.num_successors[u], self.st.num_successors[v]);
            if u_num <= v_num {
                u = self.st.nodes[u].parent;
            }
            if v_num <= u_num {
                v = self.st.nodes[v].parent;
            }
        }
        u
    }
}
