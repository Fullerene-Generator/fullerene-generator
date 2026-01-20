#include <expansions/l_reduction.h>

#include <cstdint>
#include <vector>
#include <queue>

namespace {

	using Pmask = std::uint64_t;

	static inline Pmask pbit(std::size_t i) {
		return Pmask{ 1 } << i;
	}

	static Pmask pentagons_within_4(const std::shared_ptr<base_node>& start_node,
		int total_nodes)
	{
		std::vector<std::int8_t> dist(static_cast<std::size_t>(total_nodes), -1);
		std::queue<std::shared_ptr<base_node>> q;

		const int start_id = static_cast<int>(start_node->id());
		dist[static_cast<std::size_t>(start_id)] = 0;
		q.push(start_node);

		Pmask mask = 0;

		while (!q.empty()) {
			const auto node = q.front();
			q.pop();

			const int v = static_cast<int>(node->id());
			const std::int8_t dv = dist[static_cast<std::size_t>(v)];
			if (dv > 4) continue;

			mask |= pbit(static_cast<std::size_t>(v));
			if (dv == 4) continue;

			for (const auto& wweak : node->neighbors()) {
				if (auto w = wweak.lock()) {
					const int to = static_cast<int>(w->id());
					auto& dref = dist[static_cast<std::size_t>(to)];
					if (dref == -1) {
						dref = static_cast<std::int8_t>(dv + 1);
						q.push(w);
					}
				}
			}
		}

		return mask;
	}

}


bool has_L0_pair_pent_distance_gt4(const dual_fullerene& G)
{
	const auto l0s = find_l_reductions(G, 1);
	if (l0s.size() < 2) return false;

	const auto& pent_nodes = G.get_nodes_5();
	const std::size_t P = pent_nodes.size();
	if (P < 2) return false;

	const int n = static_cast<int>(G.total_nodes());


	// near[i] = bitmask of pentagons within distance <= 4 from pentagon i
	std::vector<Pmask> near(P, 0);
	for (std::size_t i = 0; i < P; ++i) {
		near[i] = pentagons_within_4(pent_nodes[i], n);
	}

	// For each reduction store union of pentagons near the reduction pentagons
	std::vector<Pmask> prev_union;
	prev_union.reserve(l0s.size());

	for (const auto& r : l0s) {
		const int a = static_cast<int>(r.first_edge.from->id());
		const int b = static_cast<int>(r.second_edge.from->id());

		const Pmask bits_ab = pbit(static_cast<std::size_t>(a)) | pbit(static_cast<std::size_t>(b));
		for (const Pmask prev : prev_union) {
			if ((prev & bits_ab) == 0) return true;
		}

		prev_union.push_back(near[static_cast<std::size_t>(a)] | near[static_cast<std::size_t>(b)]);
	}

	return false;
}




std::vector<l_reduction>
find_l_reductions(const dual_fullerene& G, int size)
{
	std::vector<l_reduction> out;
	if (size < 1) return out;

	int path_len = size + 1;

	for (const auto& pent : G.get_nodes_5()) {
		const auto& start_node = pent;
		int deg = static_cast<int>(start_node->degree());

		for (int i = 0; i < deg; ++i) {
			directed_edge e0{ start_node, static_cast<std::size_t>(i) };

			for (bool use_next : { true, false }) {
				int outside_hex1 = use_next
					? static_cast<int>(e0.prev_around(2).to()->degree())
					: static_cast<int>(e0.next_around(2).to()->degree());
				if (outside_hex1 != 6) continue;

				auto e = e0;
				std::vector<int> path;
				path.reserve(path_len);

				bool ok = true;
				for (int step = 0; step < path_len; ++step) {
					auto v = e.from;

					if (step == 0 || step == path_len - 1) {
						if (v->type() != node_type::NODE_5) { ok = false; break; }
					}
					else {
						if (v->type() != node_type::NODE_6) { ok = false; break; }
					}

					path.push_back(static_cast<int>(v->id()));

					if (step < path_len - 1) {
						e = use_next ? e.right_turn(3) : e.left_turn(3);
					}
				}
				if (!ok) continue;

				int outside_hex2 = use_next
					? static_cast<int>(e.next_around(1).to()->degree())
					: static_cast<int>(e.prev_around(1).to()->degree());
				if (outside_hex2 != 6) continue;

				int last_pent = path.back();

				auto last_node = G.get_node(static_cast<unsigned int>(last_pent));
				auto prev_node = G.get_node(static_cast<unsigned int>(path[path.size() - 2]));
				directed_edge second_edge = last_node->get_edge(prev_node);

				l_reduction r;
				r.first_edge = e0;
				r.second_edge = second_edge;
				r.use_next = use_next;
				r.size = size;

				out.push_back(r);
			}
		}
	}

	return out;
}

void l_reduction::apply(dual_fullerene& G, const expansion_candidate& c) const
{
	const int i = size - 1;

	const int created = i + 2;
	const int h1 = static_cast<int>(G.total_nodes()) - created;
	const int h2 = static_cast<int>(G.total_nodes()) - 1;

	auto h2_node = G.get_node(static_cast<unsigned int>(h2));

	int u_first = c.path[i + 1];
	int u_second = c.path[i + 2];
	int w_first = c.parallel_path[i + 1];
	int w_second = c.parallel_path[i + 2];
	auto u_first_node = G.get_node(u_first);
	auto u_second_node = G.get_node(u_second);
	auto w_second_node = G.get_node(w_second);

	G.replace_neighbor(w_first, w_second, u_second);
	G.replace_neighbor(u_second, w_second, w_first);
	G.replace_neighbor(u_first, w_second, w_first);

	h2_node->remove_neighbor(w_second_node);
	G.move_neighborhood(h2, w_second);
	int h = h2 - 1;
	G.pop_last_node6();

	for (int j = i; j > 0; --j) {
		auto h_node = G.get_node(h);

		u_first = c.path[j];
		u_second = c.path[j + 1];
		w_first = c.parallel_path[j];
		w_second = c.parallel_path[j + 1];

		auto u_first_node2 = G.get_node(u_first);
		auto u_second_node2 = G.get_node(u_second);
		auto w_first_node2 = G.get_node(w_first);
		auto w_second_node2 = G.get_node(w_second);

		w_second_node2->replace_neighbor(h_node, u_second_node2);
		w_first_node2->replace_neighbor(h_node, u_second_node2);
		u_second_node2->replace_neighbor(h_node, w_first_node2);
		u_first_node2->replace_neighbor(h_node, w_first_node2);

		G.pop_last_node6();
		--h;
	}

	auto h1_node = G.get_node(h1);

	u_first = c.path[0];
	u_second = c.path[1];
	w_first = c.parallel_path[0];
	w_second = c.parallel_path[1];

	auto u_first_node3 = G.get_node(u_first);
	auto u_second_node3 = G.get_node(u_second);
	auto w_first_node3 = G.get_node(w_first);
	auto w_second_node3 = G.get_node(w_second);

	w_second_node3->replace_neighbor(u_first_node3, u_second_node3);
	w_first_node3->replace_neighbor(u_first_node3, u_second_node3);
	u_second_node3->replace_neighbor(u_first_node3, w_first_node3);
	h1_node->remove_neighbor(u_first_node3);

	G.move_neighborhood(h1, u_first);
	G.pop_last_node6();
}
