use network_algorithms::minimum_cost_flow::graph::Graph;
use network_algorithms::minimum_cost_flow::network_simplex_pivot_rules::{BlockSearchPivotRule, PivotRule};
use network_algorithms::minimum_cost_flow::primal_network_simplex::PrimalNetworkSimplex;
use network_algorithms::minimum_cost_flow::status::Status;

fn main() {
    let mut graph = Graph::default();
    graph.add_nodes(4);

    let edges = vec![
        graph.add_directed_edge(0, 1, 0, 2, 1).unwrap(),
        graph.add_directed_edge(0, 2, 0, 1, 2).unwrap(),
        graph.add_directed_edge(1, 2, 0, 1, 1).unwrap(),
        graph.add_directed_edge(1, 3, 0, 1, 3).unwrap(),
        graph.add_directed_edge(2, 3, 0, 2, 1).unwrap(),
    ];

    graph.add_supply(0, 2);
    graph.add_supply(3, -2);

    let status = PrimalNetworkSimplex::default().solve(&mut BlockSearchPivotRule::new(edges.len()), &mut graph);

    match status {
        Status::Optimal => {
            println!("minimum cost:{}", graph.minimum_cost());
            for edge_id in edges {
                println!("{:?}", graph.get_edge(edge_id).unwrap());
            }
        }
        _ => unreachable!(),
    }
}
