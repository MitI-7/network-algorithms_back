use crate::minimum_cost_flow::spanning_tree_structure::{EdgeState, InternalEdge, SpanningTreeStructure};
use num_traits::NumAssign;
use std::ops::Neg;

pub trait PivotRule<Flow> {
    fn new(num_edges: usize) -> Self;
    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize>;
}

pub struct BestEligibleArcPivotRule<Flow> {
    _maker: std::marker::PhantomData<fn() -> Flow>,
}

impl<Flow> PivotRule<Flow> for BestEligibleArcPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    fn new(_num_edges: usize) -> Self {
        Self { _maker: std::marker::PhantomData }
    }

    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize> {
        let mut maxi_violation = Flow::zero();
        let mut entering_edge_id = None;

        for (edge_id, edge) in st.edges.iter().enumerate() {
            let violation = calculate_violation(edge, st);
            if violation > maxi_violation {
                maxi_violation = violation;
                entering_edge_id = Some(edge_id);
            }
        }

        entering_edge_id
    }
}

pub struct FirstEligibleArcPivotRule<Flow> {
    current_edge_id: usize,
    _maker: std::marker::PhantomData<fn() -> Flow>,
}

impl<Flow> PivotRule<Flow> for FirstEligibleArcPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    fn new(_num_edges: usize) -> Self {
        Self { current_edge_id: 0, _maker: std::marker::PhantomData }
    }

    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize> {
        for _ in 0..st.num_edges {
            let edge = &st.edges[self.current_edge_id];
            let violation = calculate_violation(edge, st);

            if violation > Flow::zero() {
                return Some(self.current_edge_id);
            }

            self.current_edge_id += 1;
            if self.current_edge_id == st.num_edges {
                self.current_edge_id = 0;
            }
        }

        None
    }
}

pub struct BlockSearchPivotRule<Flow> {
    current_edge_id: usize,
    block_size: usize,
    _maker: std::marker::PhantomData<fn() -> Flow>,
}

impl<Flow> BlockSearchPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn new_with_parameter(num_edges: usize, min_block_size: usize, block_size_factor: f64) -> Self {
        assert!(min_block_size > 0);
        assert!(block_size_factor >= 0.0);
        Self { current_edge_id: 0, block_size: min_block_size.max((block_size_factor * (num_edges as f64).sqrt()) as usize), _maker: std::marker::PhantomData }
    }
}

impl<Flow> PivotRule<Flow> for BlockSearchPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    fn new(num_edges: usize) -> Self {
        let min_block_size = 10;
        let block_size_factor = 1.0; // between 0.5 and 2.0
        Self::new_with_parameter(num_edges, min_block_size, block_size_factor)
    }

    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize> {
        let mut maxi_violation = Flow::zero();
        let mut entering_edge_id = None;
        let mut count = self.block_size;

        for _ in 0..st.num_edges {
            let edge = &st.edges[self.current_edge_id];
            let violation = calculate_violation(edge, st);

            if violation > maxi_violation {
                maxi_violation = violation;
                entering_edge_id = Some(self.current_edge_id);
            }

            count -= 1;
            if count == 0 {
                if entering_edge_id.is_some() {
                    return entering_edge_id;
                }
                count = self.block_size;
            }

            self.current_edge_id += 1;
            if self.current_edge_id == st.num_edges {
                self.current_edge_id = 0;
            }
        }

        entering_edge_id
    }
}

pub struct CandidateListPivotRule<Flow> {
    current_edge_id: usize,
    candidates: Box<[usize]>,
    candidate_list_size: usize,
    minor_count_limit: usize,
    minor_count: usize,
    current_size: usize,
    _maker: std::marker::PhantomData<fn() -> Flow>,
}

impl<Flow> CandidateListPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn new_with_parameter(num_edges: usize, min_candidate_list_size: usize, candidate_list_size_factor: f64, min_minor_limit: usize, minor_limit_factor: f64) -> Self {
        assert!(min_candidate_list_size > 0);
        assert!(candidate_list_size_factor > 0.0);
        assert!(min_minor_limit > 0);
        assert!(minor_limit_factor >= 0.0);

        let candidate_list_size = min_candidate_list_size.max((candidate_list_size_factor * (num_edges as f64).sqrt()) as usize);
        let minor_limit = min_minor_limit.max((minor_limit_factor * candidate_list_size as f64) as usize);

        Self {
            current_edge_id: 0,
            candidates: vec![usize::MAX; candidate_list_size].into_boxed_slice(),
            candidate_list_size,
            current_size: 0,
            minor_count_limit: minor_limit,
            minor_count: 0,
            _maker: std::marker::PhantomData,
        }
    }
}

impl<Flow> PivotRule<Flow> for CandidateListPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    fn new(num_edges: usize) -> Self {
        let min_candidate_list_size = 10;
        let candidate_list_size_factor = 0.25;
        let min_minor_limit = 3;
        let minor_limit_factor = 0.1;
        Self::new_with_parameter(num_edges, min_candidate_list_size, candidate_list_size_factor, min_minor_limit, minor_limit_factor)
    }

    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize> {
        let mut maxi_violation = Flow::zero();
        let mut entering_edge_id = None;

        // minor iteration
        if self.current_size > 0 && self.minor_count < self.minor_count_limit {
            self.minor_count += 1;

            // search in candidate list
            let mut i = 0;
            while i < self.current_size {
                let edge_id = self.candidates[i];
                let edge = &st.edges[edge_id];
                let violation = calculate_violation(edge, st);

                if violation <= Flow::zero() {
                    // remove ineligible arc from the candidates
                    self.current_size -= 1;
                    self.candidates[i] = self.candidates[self.current_size];
                } else {
                    if violation > maxi_violation {
                        maxi_violation = violation;
                        entering_edge_id = Some(edge_id);
                    }
                    i += 1;
                }
            }

            if entering_edge_id.is_some() {
                return entering_edge_id;
            }
        }

        // build a candidate list
        self.current_size = 0;
        for _ in 0..st.num_edges {
            let edge = &st.edges[self.current_edge_id];
            let violation = match edge.state {
                EdgeState::Upper => st.reduced_cost(edge),
                _ => -st.reduced_cost(edge),
            };

            if violation > Flow::zero() {
                self.candidates[self.current_size] = self.current_size;
                self.current_size += 1;

                if violation > maxi_violation {
                    maxi_violation = violation;
                    entering_edge_id = Some(self.current_edge_id);
                }
            }

            if self.current_size == self.candidate_list_size {
                break;
            }

            self.current_edge_id += 1;
            if self.current_edge_id == st.num_edges {
                self.current_edge_id = 0;
            }
        }

        self.minor_count = 1;
        entering_edge_id
    }
}

pub struct AlteringCandidateListPivotRule<Flow> {
    current_edge_id: usize,
    block_size: usize,
    head_length: usize,
    candidates: Box<[(usize, Flow)]>,
    current_size: usize,
}

impl<Flow> AlteringCandidateListPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    pub fn new_with_parameter(num_edges: usize, min_block_size: usize, block_size_factor: f64, min_head_length: usize, head_length_factor: f64) -> Self {
        assert!(min_block_size > 0);
        assert!(block_size_factor > 0.0);
        assert!(min_head_length > 0);
        assert!(head_length_factor >= 0.0);

        let block_size = min_block_size.max((block_size_factor * (num_edges as f64).sqrt()) as usize);
        let head_length = min_head_length.max((head_length_factor * block_size as f64) as usize);

        Self { current_edge_id: 0, block_size, head_length, candidates: vec![(usize::MAX, Flow::zero()); head_length + block_size].into_boxed_slice(), current_size: 0 }
    }
}

impl<Flow> PivotRule<Flow> for AlteringCandidateListPivotRule<Flow>
where
    Flow: NumAssign + Neg<Output = Flow> + Ord + Copy,
{
    fn new(num_edges: usize) -> Self {
        let min_block_size = 10;
        let block_size_factor = 1.0;
        let min_head_length = 3;
        let head_length_factor = 0.01;
        Self::new_with_parameter(num_edges, min_block_size, block_size_factor, min_head_length, head_length_factor)
    }

    fn find_entering_edge<F: Fn(&InternalEdge<Flow>, &SpanningTreeStructure<Flow>) -> Flow>(&mut self, st: &SpanningTreeStructure<Flow>, calculate_violation: F) -> Option<usize> {
        // update candidate cost
        let mut i = 0;
        while i < self.current_size {
            let (edge_id, _) = self.candidates[i];
            let edge = &st.edges[edge_id];
            let violation = calculate_violation(edge, st);

            if violation <= Flow::zero() {
                // remove ineligible arc from the candidates
                self.current_size -= 1;
                self.candidates[i] = self.candidates[self.current_size];
            } else {
                self.candidates[i].1 = violation;
                i += 1;
            }
        }

        // extend the candidate list
        let mut block_count = self.block_size;
        let mut limit = self.head_length;

        for _ in 0..st.num_edges {
            let edge = &st.edges[self.current_edge_id];
            let violation = match edge.state {
                EdgeState::Upper => st.reduced_cost(edge),
                _ => -st.reduced_cost(edge),
            };

            // add eligible arc to the candidates
            if violation > Flow::zero() {
                self.candidates[self.current_size] = (self.current_edge_id, violation);
                self.current_size += 1;
            }
            block_count -= 1;

            if block_count == 0 {
                if self.current_size > limit {
                    break;
                }
                limit = 0;
                block_count = self.block_size;
            }

            self.current_edge_id += 1;
            if self.current_edge_id == st.num_edges {
                self.current_edge_id = 0;
            }
        }

        if self.current_size == 0 {
            return None;
        }

        let new_length = self.current_size.min(self.head_length + 1);
        if new_length == self.current_size {
            self.candidates[..self.current_size].sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
        } else {
            self.candidates[..self.current_size].select_nth_unstable_by(new_length, |a, b| b.1.partial_cmp(&a.1).unwrap());
        }

        let entering_edge_id = Some(self.candidates[0].0);
        self.candidates[0] = self.candidates[new_length - 1];
        self.current_size = new_length - 1;

        entering_edge_id
    }
}
