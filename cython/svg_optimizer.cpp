#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <cstdio>

#define INF 2000000000

using Level = int32_t;
using Index = int32_t;
using ID = int32_t;
using Position = int32_t;

const std::vector<std::vector<Position>> permutations = {
    {0, 1, 2, 3}, {3, 0, 1, 2}, {2, 3, 0, 1}, {1, 2, 3, 0},
    {2, 1, 0, 3}, {3, 2, 1, 0}, {0, 3, 2, 1}, {1, 0, 3, 2}
};

struct Solution {
    std::vector<Position> positions;
    int score = -1;
    std::vector<std::pair<ID, ID>> dangling_edges;
};

Level LevelOf(ID id) { return id / 4; }

int Score(const std::pair<Level, Position> &a,
          const std::pair<Level, Position> &b) {

    if (a.first == b.first) return 0;

    if (a.second == b.second && std::abs(a.first - b.first) == 1)
        return 0;

    if ((a.second == 2 && b.second == 2) || (a.second == 3 && b.second == 3))
        return 3;

    if ((a.second == 2 && b.second == 3) || (a.second == 3 && b.second == 2))
        return 5;

    if ((a.second == 0 && b.second == 0) || (a.second == 1 && b.second == 1))
        return 3;

    if ((a.second == 0 && b.second == 1) || (a.second == 1 && b.second == 0))
        return 5;

    if ((a.second == 0 && b.second == 3) || (a.second == 3 && b.second == 0))
        return 7;

    if ((a.second == 1 && b.second == 2) || (a.second == 2 && b.second == 1))
        return 4;

    if ((a.second == 1 && b.second == 3) ||
        (a.second == 3 && b.second == 1) ||
        (a.second == 0 && b.second == 2) ||
        (a.second == 2 && b.second == 0))
        return 500;

    return 500;
}

void UpdateScoreForLevel(const std::vector<ID> &edges, Level level,
                         Solution *solution) {

    auto score = [&](ID start, ID end) {
        return Score({LevelOf(start), solution->positions[start]},
                     {LevelOf(end), solution->positions[end]});
    };

    for (auto &[start, end] : solution->dangling_edges) {
        if (LevelOf(end) == level)
            solution->score += score(start, end);
    }

    auto &d = solution->dangling_edges;
    std::erase_if(d, [level](const auto &e) {
        return LevelOf(e.second) == level;
    });

    for (Index i = 0; i < 4; ++i) {
        ID start = level * 4 + i;
        ID end = edges[start];

        if (LevelOf(end) > level)
            solution->dangling_edges.emplace_back(start, end);
        else if (end != -1)
            solution->score += score(start, end);
    }
}

/* =========================
   [ZMIANA] FIX SameLayout
   ========================= */
bool SameLayout(const std::vector<Position> &previous,
                const std::vector<Position> &current,
                const std::vector<ID> &alignments,
                Level level) {

    for (size_t i = 0; i < 4; i++) {

        if (alignments[previous[i] + (level - 1) * 4] !=
            alignments[current[i] + level * 4]) {
            return false;
        }
    }
    return true;
}

/* =========================
   [ZMIANA] Solve signature FIX
   ========================= */
Solution Solve(const std::vector<ID> &edges,
               const std::vector<Level> &rotations,
               const std::vector<ID> &alignments) {

    const int num_levels = edges.size() / 4;
    std::vector<Solution> current, next;

    for (const auto &perm : permutations) {
        current.push_back({perm, 0});
        UpdateScoreForLevel(edges, 0, &current.back());
    }

    int best_score = INF;

    for (int level = 1; level < num_levels; ++level) {

        next.clear();
        best_score = INF;

        for (const auto &prev_sol : current) {
            for (const auto &perm : permutations) {

                if (rotations[level] != -1 &&
                    rotations[level] == rotations[level - 1]) {

                    std::vector<Position> prev(prev_sol.positions.end() - 4,
                                                prev_sol.positions.end());

                    if (!SameLayout(prev, perm, alignments, level))
                        continue;
                }

                auto &next_sol = next.emplace_back(prev_sol);

                next_sol.positions.insert(next_sol.positions.end(),
                                          perm.begin(), perm.end());

                UpdateScoreForLevel(edges, level, &next_sol);

                if (next_sol.score < best_score)
                    best_score = next_sol.score;
            }
        }

        if (!next.empty()) {
            int worst_case = best_score + next.back().dangling_edges.size() * 20;

            std::erase_if(next, [worst_case](const auto &e) {
                return e.score > worst_case;
            });
        }

        std::swap(current, next);
    }

    Solution best;
    best.score = INF;

    for (const auto &s : current)
        if (s.score < best.score)
            best = s;

    return best;
}

/* =========================
   [ZMIANA] FIX SolveFailsafe
   ========================= */
Solution SolveFailsafe(const std::vector<ID> &edges,
                       const std::vector<Level> &rotations,
                       const std::vector<ID> &alignments) {

    Solution result = Solve(edges, rotations, alignments);

    if (result.score == -1) {
        std::vector<Level> rot_fixed(edges.size(), -1);
        result = Solve(edges, rot_fixed, alignments);
    }

    return result;
}

/* =========================
   [ZMIANA] FIX main indexing
   ========================= */
int main() {

    std::vector<ID> edges;
    std::vector<Level> rotations;
    std::vector<ID> alignments;

    int nucl_num, tmp;

    scanf("%d", &nucl_num);

    for (int i = 0; i < nucl_num; ++i) {
        scanf("%d", &tmp);
        edges.push_back(tmp);
    }

    for (int i = 0; i < nucl_num / 4; ++i) {
        scanf("%d", &tmp);
        rotations.push_back(tmp);
    }

    for (int i = 0; i < nucl_num; ++i) {
        scanf("%d", &tmp);
        alignments.push_back(tmp);
    }

    Solution result = Solve(edges, rotations, alignments);

    if (result.score == -1) {
        for (auto &r : rotations) r = -1;
        result = Solve(edges, rotations, alignments);
    }

    printf("%d", result.positions[0]);
    for (size_t i = 1; i < result.positions.size(); ++i)
        printf(" %d", result.positions[i]);

    return 0;
}
