use crate::minimum_cost_flow::graph::Graph;
use num_traits::NumAssign;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::ops::Neg;

#[derive(Default, Clone)]
pub struct Node<Flow> {
    pub parent: usize,
    pub parent_edge_id: usize,
    pub potential: Flow,
}

#[derive(Default, PartialEq, Debug)]
pub enum EdgeState {
    #[default]
    Lower,
    Upper,
    Tree,
}

#[derive(Default)]
pub struct InternalEdge<Flow> {
    pub from: usize,
    pub to: usize,
    pub upper: Flow,
    pub cost: Flow,
    pub flow: Flow,
    pub state: EdgeState,
}

impl<Flow> InternalEdge<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn is_feasible(&self) -> bool {
        Flow::zero() <= self.flow && self.flow <= self.upper
    }

    pub fn is_lower(&self) -> bool {
        self.flow == Flow::zero()
    }

    pub fn is_upper(&self) -> bool {
        self.flow == self.upper
    }

    pub fn residual_capacity(&self) -> Flow {
        self.upper - self.flow
    }

    pub fn opposite_side(&self, u: usize) -> usize {
        debug_assert!(u == self.from || u == self.to);
        u ^ self.to ^ self.from
    }
}

#[derive(Default)]
pub struct SpanningTreeStructure<Flow> {
    pub(crate) num_nodes: usize,
    pub(crate) num_edges: usize,
    pub(crate) excesses: Vec<Flow>,

    pub(crate) nodes: Vec<Node<Flow>>,
    pub(crate) edges: Vec<InternalEdge<Flow>>,

    pub(crate) root: usize,
    pub(crate) next_node_dft: Vec<usize>,       // next nodes in depth-first thread
    pub(crate) prev_node_dft: Vec<usize>,       // previous nodes in depth-first thread
    pub(crate) last_descendent_dft: Vec<usize>, // last descendants in depth-first thread
    pub(crate) num_successors: Vec<usize>,      // the number of successors of the node in the tree
}

#[allow(dead_code)]
impl<Flow> SpanningTreeStructure<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy + Clone,
{
    pub(crate) fn build(&mut self, graph: &mut Graph<Flow>) {
        (self.num_nodes, self.num_edges) = (graph.num_nodes(), graph.num_edges());
        self.excesses = graph.excesses.clone();

        for edge in graph.edges.iter() {
            assert!(edge.upper >= Flow::zero() && edge.cost >= Flow::zero());
            self.edges
                .push(InternalEdge { from: edge.from, to: edge.to, flow: edge.flow, upper: edge.upper, cost: edge.cost, state: EdgeState::Lower });
        }

        self.root = usize::MAX;
        self.nodes.resize(self.num_nodes, Node { parent: usize::MAX, parent_edge_id: usize::MAX, potential: Flow::zero() });
        self.next_node_dft.resize(self.num_nodes, usize::MAX);
        self.prev_node_dft.resize(self.num_nodes, usize::MAX);
        self.last_descendent_dft.resize(self.num_nodes, usize::MAX);
        self.num_successors.resize(self.num_nodes, 0);
    }

    #[inline]
    pub(crate) fn reduced_cost(&self, edge: &InternalEdge<Flow>) -> Flow {
        edge.cost - self.nodes[edge.from].potential + self.nodes[edge.to].potential
    }

    pub(crate) fn update_flow_in_path(&mut self, source: usize, sink: usize, delta: Flow) {
        let mut now = sink;
        while now != source {
            let (parent, edge_id) = (self.nodes[now].parent, self.nodes[now].parent_edge_id);
            let edge = &mut self.edges[edge_id];
            edge.flow += if edge.from == parent { delta } else { -delta };
            now = parent;
        }
        self.excesses[source] -= delta;
        self.excesses[sink] += delta;
    }

    pub(crate) fn update_flow_in_cycle(&mut self, entering_edge_id: usize, delta: Flow, apex: usize) {
        let delta = match self.edges[entering_edge_id].state {
            EdgeState::Upper => -delta,
            _ => delta,
        };
        self.edges[entering_edge_id].flow += delta;

        let mut now = self.edges[entering_edge_id].from;
        while now != apex {
            let edge = &mut self.edges[self.nodes[now].parent_edge_id];
            edge.flow += if now == edge.from { -delta } else { delta };
            now = self.nodes[now].parent;
        }

        let mut now = self.edges[entering_edge_id].to;
        while now != apex {
            let edge = &mut self.edges[self.nodes[now].parent_edge_id];
            edge.flow += if now == edge.from { delta } else { -delta };
            now = self.nodes[now].parent;
        }
    }

    // change the root of subtree from now_root to new_root
    // O(|tree|)
    pub(crate) fn re_rooting(&mut self, _now_root: usize, new_root: usize, entering_edge_id: usize) {
        let mut ancestors = Vec::new();
        let mut now = new_root;
        while now != usize::MAX {
            ancestors.push(now);
            now = self.nodes[now].parent;
        }
        ancestors.reverse();

        for pair in ancestors.windows(2) {
            let (p, q) = (pair[0], pair[1]);
            let size_p = self.num_successors[p];
            let last_q = self.last_descendent_dft[q];

            self.nodes[p].parent = q;
            self.nodes[q].parent = usize::MAX;
            self.nodes[p].parent_edge_id = self.nodes[q].parent_edge_id;
            self.nodes[q].parent_edge_id = usize::MAX;
            self.num_successors[p] = size_p - self.num_successors[q];
            self.num_successors[q] = size_p;

            let prev_q = self.prev_node_dft[q];
            let next_last_q = self.next_node_dft[last_q];
            self.next_node_dft[prev_q] = next_last_q;
            self.prev_node_dft[next_last_q] = prev_q;
            self.next_node_dft[last_q] = q;
            self.prev_node_dft[q] = last_q;

            let mut last_p = self.last_descendent_dft[p];
            if last_p == last_q {
                self.last_descendent_dft[p] = prev_q;
                last_p = prev_q;
            }

            self.prev_node_dft[p] = last_q;
            self.next_node_dft[last_q] = p;
            self.next_node_dft[last_p] = q;
            self.prev_node_dft[q] = last_p;
            self.last_descendent_dft[q] = last_p;
        }

        // update potential
        let entering_edge = &self.edges[entering_edge_id];
        let delta = if new_root == entering_edge.from {
            self.reduced_cost(entering_edge)
        } else {
            -self.reduced_cost(entering_edge)
        };

        let mut now = new_root;
        while now != usize::MAX {
            self.nodes[now].potential += delta;
            if now == self.last_descendent_dft[new_root] {
                break;
            }
            now = self.next_node_dft[now];
        }
    }

    // remove leaving_edge_id
    pub(crate) fn detach_tree(&mut self, _root: usize, sub_tree_root: usize, leaving_edge_id: usize) {
        let leaving_edge = &mut self.edges[leaving_edge_id];
        leaving_edge.state = if leaving_edge.is_lower() { EdgeState::Lower } else { EdgeState::Upper };

        // detach sub tree
        self.nodes[sub_tree_root].parent = usize::MAX;
        self.nodes[sub_tree_root].parent_edge_id = usize::MAX;

        let prev_t = self.prev_node_dft[sub_tree_root];
        let last_t = self.last_descendent_dft[sub_tree_root];
        let next_last_t = self.next_node_dft[last_t];
        self.next_node_dft[prev_t] = next_last_t;
        self.prev_node_dft[next_last_t] = prev_t;
        self.next_node_dft[last_t] = sub_tree_root;
        self.prev_node_dft[sub_tree_root] = last_t;

        let sub_tree_size = self.num_successors[sub_tree_root];
        let mut now = leaving_edge.opposite_side(sub_tree_root);
        while now != usize::MAX {
            self.num_successors[now] -= sub_tree_size;
            if self.last_descendent_dft[now] == last_t {
                self.last_descendent_dft[now] = prev_t;
            }
            now = self.nodes[now].parent;
        }
    }

    // attach T2 under T1
    // O(1)
    // add entering_ege_id
    pub(crate) fn attach_tree(&mut self, _root: usize, attach_node: usize, sub_tree_root: usize, entering_edge_id: usize) {
        self.edges[entering_edge_id].state = EdgeState::Tree;

        let (p, q) = (attach_node, sub_tree_root); // p -> q

        // attach tree
        self.nodes[q].parent = p;
        self.nodes[q].parent_edge_id = entering_edge_id;

        let last_p = self.last_descendent_dft[attach_node];
        let next_last_p = self.next_node_dft[last_p];
        let last_q = self.last_descendent_dft[q];
        self.next_node_dft[last_p] = q;
        self.prev_node_dft[q] = last_p;
        self.prev_node_dft[next_last_p] = last_q;
        self.next_node_dft[last_q] = next_last_p;

        let sub_tree_size = self.num_successors[q];
        let mut now = attach_node;
        while now != usize::MAX {
            self.num_successors[now] += sub_tree_size;
            if self.last_descendent_dft[now] == last_p {
                self.last_descendent_dft[now] = last_q
            }
            now = self.nodes[now].parent;
        }
    }

    // dijkstra
    pub(crate) fn shortest_path(&self, source: usize) -> (Vec<Flow>, Vec<Option<usize>>) {
        let mut graph = vec![Vec::new(); self.num_nodes];
        let mut total_cost = Flow::zero();
        for (edge_id, edge) in self.edges.iter().enumerate() {
            graph[edge.from].push(edge_id);
            assert!(edge.cost >= Flow::zero());
            total_cost += edge.cost;
        }

        let mut distances = vec![total_cost + Flow::one(); self.num_nodes];
        let mut prev_edge_id = vec![None; self.num_nodes];
        let mut seen = vec![false; self.num_nodes];
        let mut bh = BinaryHeap::from([(Reverse(Flow::zero()), source)]);

        distances[source] = Flow::zero();
        while let Some((now_dist, u)) = bh.pop() {
            if seen[u] {
                continue;
            }
            seen[u] = true;

            for &edge_id in graph[u].iter() {
                let edge = &self.edges[edge_id];
                let new_dist = now_dist.0 + edge.cost;

                if new_dist < distances[edge.to] {
                    prev_edge_id[edge.to] = Some(edge_id);
                    distances[edge.to] = new_dist;
                    bh.push((Reverse(new_dist), edge.to));
                }
            }
        }

        (distances, prev_edge_id)
    }

    pub fn satisfy_constraints(&self) -> bool {
        self.edges.iter().all(|edge| edge.is_feasible()) && self.excesses.iter().all(|&excess| excess == Flow::zero())
    }

    pub fn satisfy_optimality_conditions(&self) -> bool {
        self.edges.iter().all(|edge| match edge.state {
            EdgeState::Tree => self.reduced_cost(edge) == Flow::zero(),
            EdgeState::Lower => edge.upper == Flow::zero() || self.reduced_cost(edge) >= Flow::zero(),
            EdgeState::Upper => edge.upper == Flow::zero() || self.reduced_cost(edge) <= Flow::zero(),
        })
    }

    pub fn validate_num_successors(&self, root: usize) -> bool {
        let mut order = Vec::new();
        let mut now = root;
        loop {
            order.push(now);
            now = self.next_node_dft[now];
            if now == root {
                break;
            }
        }

        let mut num_successors = vec![1; self.num_nodes];
        for &u in order.iter().rev() {
            if num_successors[u] != self.num_successors[u] {
                return false;
            }
            if self.nodes[u].parent != usize::MAX {
                num_successors[self.nodes[u].parent] += num_successors[u];
            }
        }

        true
    }
}
