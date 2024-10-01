#[derive(Default, PartialEq, Debug)]
pub enum Status {
    #[default]
    NotSolved,
    BadInput,
    Unbalanced,
    Infeasible,
    Optimal,
}
