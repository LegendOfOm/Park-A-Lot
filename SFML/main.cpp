#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <optional>
#include <algorithm>

using std::cout;
using std::cin;
using std::vector;

// --------------------------------------------------------
// Geometry / lot constants
// --------------------------------------------------------

constexpr double STALL_WIDTH  = 2.5;  // along lane (frontage)
constexpr double STALL_DEPTH  = 5.0;  // perpendicular to lane

// entrance/exit lane width = stall width + 0.5m, as you requested
constexpr double LANE_W       = STALL_WIDTH + 0.5;

// --------------------------------------------------------
// Basic types
// --------------------------------------------------------

struct Rect {
    double x0, y0;
    double x1, y1; // x1 > x0, y1 > y0
};

struct Stall {
    Rect r;
};

struct Lane {
    Rect r;
};

enum class Edge { Bottom, Top, Left, Right };

enum class StallPattern {
    Perpendicular, // classic 90° parking (what you have now)
    Parallel       // long side along the lane (parallel parking)
};


struct Gate {
    Edge   edge;     // which side of lot
    double offset;   // along that edge (metres)
    double width;    // opening width (metres)
};

struct Layout {
    bool valid = false;

    vector<Stall> stalls;   // all accessible stalls
    vector<Lane>  lanes;    // continuous snake path

    Gate entrance{};
    Gate exit{};

    int    stallCount   = 0;
    double usedArea     = 0.0;
    double wastedArea   = 0.0;
};

struct ParkingRegion {
    Rect area;   // free rectangle
    Edge side;   // which side is touched by a lane (Bottom, Top, Left, Right)
};


// --------------------------------------------------------
// Rect helpers
// --------------------------------------------------------

double rectArea(const Rect& r)
{
    return (r.x1 - r.x0) * (r.y1 - r.y0);
}

// --------------------------------------------------------
// Vertical snake generator (lanes run along Y)
// Stalls are perpendicular to lanes.
// --------------------------------------------------------

// --------------------------------------------------------
// Overlap helpers
// --------------------------------------------------------

bool rectsOverlap(const Rect& a, const Rect& b)
{
    double xOverlap = std::min(a.x1, b.x1) - std::max(a.x0, b.x0);
    double yOverlap = std::min(a.y1, b.y1) - std::max(a.y0, b.y0);
    return xOverlap > 1e-6 && yOverlap > 1e-6;
}

// Clean up layout: remove stalls that overlap lanes and recompute areas.
void finalizeLayout(Layout& L, double lotW, double lotL)
{
    if (!L.valid) return;

    const double lotArea = lotW * lotL;

    // Lane area is always fully used
    double laneArea = 0.0;
    for (const auto& ln : L.lanes)
        laneArea += rectArea(ln.r);

    // Keep only stalls that don't collide with any lane
    std::vector<Stall> cleaned;
    cleaned.reserve(L.stalls.size());
    double stallArea = 0.0;

    for (const auto& s : L.stalls)
    {
        bool bad = false;
        for (const auto& ln : L.lanes)
        {
            if (rectsOverlap(s.r, ln.r))
            {
                bad = true;
                break;
            }
        }
        if (!bad)
        {
            cleaned.push_back(s);
            stallArea += rectArea(s.r);
        }
    }

    L.stalls     = std::move(cleaned);
    L.stallCount = static_cast<int>(L.stalls.size());
    L.usedArea   = laneArea + stallArea;
    L.wastedArea = std::max(0.0, lotArea - L.usedArea);
}

// Ensure we have two full-height vertical connectors (for short lots
// where only one vertical lane was created).
void ensureTwoVerticalConnectors(Layout& L, double lotW, double lotL)
{
    if (!L.valid || L.lanes.empty()) return;

    // Overall vertical range covered by lanes
    double minY = lotL;
    double maxY = 0.0;
    for (const auto& ln : L.lanes) {
        minY = std::min(minY, ln.r.y0);
        maxY = std::max(maxY, ln.r.y1);
    }

    // Find lanes that are "mostly vertical" and span the full height
    std::vector<int> fullVertIdx;
    for (int i = 0; i < static_cast<int>(L.lanes.size()); ++i) {
        const Rect& r = L.lanes[i].r;
        double w = r.x1 - r.x0;
        double h = r.y1 - r.y0;

        bool spansFullHeight = (r.y0 <= minY + 1e-6) && (r.y1 >= maxY - 1e-6);
        bool isVerticalStrip = (w <= LANE_W + 1e-3) && (h > w); // tall & skinny

        if (spansFullHeight && isVerticalStrip) {
            fullVertIdx.push_back(i);
        }
    }

    // If there is exactly ONE full-height vertical connector, mirror it
    // on the opposite side so we get a proper loop.
    if (fullVertIdx.size() == 1) {
        const Rect& r0 = L.lanes[fullVertIdx[0]].r;

        bool onLeft = (r0.x0 < lotW * 0.5);
        double newX0, newX1;

        if (onLeft) {
            // existing connector is left-ish -> add one on the right edge
            newX1 = lotW;
            newX0 = newX1 - LANE_W;
        } else {
            // existing connector is right-ish -> add one on the left edge
            newX0 = 0.0;
            newX1 = LANE_W;
        }

        Rect r{ newX0, minY, newX1, maxY };
        L.lanes.push_back({ r });
    }
}

bool hasTooNarrowLane(const Layout& L)
{
    for (const auto& ln : L.lanes) {
        double dx = ln.r.x1 - ln.r.x0;
        double dy = ln.r.y1 - ln.r.y0;
        double laneWidth = std::min(dx, dy);   // physical lane width

        if (laneWidth + 1e-6 < LANE_W) {
            return true;   // this lane is narrower than allowed
        }
    }
    return false;
}

std::vector<Stall> generateStallsForStrip(const ParkingRegion& R, StallPattern pattern)
{
    std::vector<Stall> out;

    double W = R.area.x1 - R.area.x0; // strip width (x)
    double H = R.area.y1 - R.area.y0; // strip height (y)

    // We assume lane runs vertically next to this strip,
    // so "along lane" is y, "away from lane" is x.

    if (pattern == StallPattern::Perpendicular) {
        // 90° parking: depth -> x, width -> y
        if (W < STALL_DEPTH || H < STALL_WIDTH) return out;

        int rowsDeep = 1; // only one row touching the lane side to be safe
        int along    = static_cast<int>(std::floor(H / STALL_WIDTH));

        double stripX0;
        if (R.side == Edge::Left) {
            // lane is on left side -> stalls start right next to lane
            stripX0 = R.area.x0;
        } else {
            // lane on right side
            stripX0 = R.area.x1 - STALL_DEPTH;
        }

        for (int i = 0; i < along; ++i) {
            double y0 = R.area.y0 + i * STALL_WIDTH;
            double y1 = y0 + STALL_WIDTH;
            double x0 = stripX0;
            double x1 = stripX0 + STALL_DEPTH;
            out.push_back({ Rect{x0, y0, x1, y1} });
        }
    }
    else if (pattern == StallPattern::Parallel) {
        // Parallel parking: width -> x (along lane), depth -> y
        if (W < STALL_WIDTH || H < STALL_DEPTH) return out;

        int along    = static_cast<int>(std::floor(W / STALL_WIDTH));
        int rowsDeep = 1; // again, one row touching the lane

        double stripY0;
        if (R.side == Edge::Left || R.side == Edge::Right) {
            // lane runs vertically, so "away" is ±y; choose side closest to lane.
            // If lane is to the left, stalls touch same y-range; here we just use bottom.
            stripY0 = R.area.y0;
        } else {
            // not used for vertical strips
            return out;
        }

        for (int i = 0; i < along; ++i) {
            double x0 = R.area.x0 + i * STALL_WIDTH;
            double x1 = x0 + STALL_WIDTH;
            double y0 = stripY0;
            double y1 = stripY0 + STALL_DEPTH;
            out.push_back({ Rect{x0, y0, x1, y1} });
        }
    }

    return out;
}

std::vector<Stall> bestStallsForRegion(const ParkingRegion& R)
{
    std::vector<Stall> best;
    int bestCount = 0;

    for (StallPattern p : {StallPattern::Perpendicular, StallPattern::Parallel}) {
        auto candidate = generateStallsForStrip(R, p);
        if (static_cast<int>(candidate.size()) > bestCount) {
            bestCount = static_cast<int>(candidate.size());
            best      = std::move(candidate);
        }
    }
    return best;
}


Layout buildVerticalSnake(double lotW, double lotL, int N)
{
    Layout L;
    L.valid = false;
    if (N <= 0) return L;

    const double moduleW = 2.0 * STALL_DEPTH + LANE_W; // [stall][lane][stall]
    const double neededW = N * moduleW;

    if (neededW > lotW + 1e-9) {
        return L; // pattern can't fit for this N
    }

    const double marginX = 0.5 * (lotW - neededW);

    // lane X positions
    vector<double> laneX0(N), laneX1(N);
    for (int i = 0; i < N; ++i) {
        laneX0[i] = marginX + STALL_DEPTH + i * moduleW;
        laneX1[i] = laneX0[i] + LANE_W;
    }

    vector<Lane> lanes;
    lanes.reserve(N + std::max(0, N - 1) + 1); // +1 for closing connector

    // vertical lanes (full height)
    for (int i = 0; i < N; ++i) {
        Rect r{ laneX0[i], 0.0, laneX1[i], lotL };
        lanes.push_back({r});
    }

    std::vector<ParkingRegion> regions;

    // left margin
    if (laneX0[0] > 0.5) {
        regions.push_back(ParkingRegion{
            Rect{0.0, 0.0, laneX0[0], lotL},
            Edge::Right  // region is touched by lane on the right
        });
    }

    // between lanes
    for (int i = 0; i < N - 1; ++i) {
        double x0 = laneX1[i];
        double x1 = laneX0[i + 1];
        if (x1 > x0 + 0.1) {
            regions.push_back(ParkingRegion{
                Rect{x0, 0.0, x1, lotL},
                Edge::Left  // or Right; you can decide a convention
            });
        }
    }

    // right margin
    if (laneX1[N-1] < lotW - 0.5) {
        regions.push_back(ParkingRegion{
            Rect{laneX1[N-1], 0.0, lotW, lotL},
            Edge::Left
        });
    }


    // ----------------------------------------------------
    // Snake connectors between adjacent lanes
    // ----------------------------------------------------
    for (int i = 0; i < N - 1; ++i) {
        double x0 = laneX1[i];
        double x1 = laneX0[i + 1];

        if (i % 2 == 0) {
            // bottom connector
            Rect r{ x0, 0.0, x1, LANE_W };
            lanes.push_back({r});
        } else {
            // top connector
            Rect r{ x0, lotL - LANE_W, x1, lotL };
            lanes.push_back({r});
        }
    }

    // ----------------------------------------------------
    // Closing connector: last lane back to first (bottom).
    // This turns the snake into a continuous loop.
    // ----------------------------------------------------
    if (N >= 2) {
        double x0 = laneX1[N - 1];
        double x1 = laneX0[0];
        if (x1 < x0) std::swap(x0, x1);
        Rect closing{ x0, 0.0, x1, LANE_W };
        lanes.push_back({closing});
    }

    // ----------------------------------------------------
    // Gates: both on bottom edge, at left and right ends.
    // They sit on the loop so you can circle indefinitely.
    // ----------------------------------------------------
    Gate entrance{};
    Gate exit{};

    entrance.edge   = Edge::Bottom;
    entrance.offset = 0.0;
    entrance.width  = LANE_W;

    exit.edge       = Edge::Bottom;
    exit.offset     = lotW - LANE_W;
    exit.width      = LANE_W;

    // ----------------------------------------------------
    // Stalls: built in stripes next to lanes, as before.
    // Every stall touches a lane, orientation 90°.
    // ----------------------------------------------------

    const int stallsPerColumn = static_cast<int>(std::floor(lotL / STALL_WIDTH));
    if (stallsPerColumn <= 0) return L; // too short to fit a single stall

    vector<Stall> stalls;
    stalls.reserve(N * stallsPerColumn * 3); // rough guess

    auto addRowOfStalls = [&](double x0, double x1) {
        for (int j = 0; j < stallsPerColumn; ++j) {
            double y0 = j * STALL_WIDTH;
            double y1 = y0 + STALL_WIDTH;
            stalls.push_back({ Rect{x0, y0, x1, y1} });
        }
    };

    // Left edge stalls (single-sided, facing lane 0)
    {
        double sx0 = marginX;
        double sx1 = laneX0[0]; // STALL_DEPTH wide
        addRowOfStalls(sx0, sx1);
    }

    // Interior double-sided strips between lanes
    for (int i = 0; i < N - 1; ++i) {
        // stripe between lane i and lane i+1 is width 2*STALL_DEPTH

        // row facing lane i (to the right)
        double rightRowX0 = laneX1[i];
        double rightRowX1 = laneX1[i] + STALL_DEPTH;
        addRowOfStalls(rightRowX0, rightRowX1);

        // row facing lane i+1 (to the left)
        double leftRowX0 = laneX0[i + 1] - STALL_DEPTH;
        double leftRowX1 = laneX0[i + 1];
        addRowOfStalls(leftRowX0, leftRowX1);
    }

    // Right edge stalls (single-sided, facing last lane)
    {
        double sx0 = laneX1[N - 1];
        double sx1 = sx0 + STALL_DEPTH;
        addRowOfStalls(sx0, sx1);
    }

    // ----------------------------------------------------
    // compute areas & fill Layout (same as before)
    // NOTE: areas still "approx" because stalls overlap
    // lanes visually. We can clean that up later.
    // ----------------------------------------------------

    const double lotArea = lotW * lotL;

    L.valid       = true;
    L.lanes       = std::move(lanes);
    L.stalls      = std::move(stalls);
    L.entrance    = entrance;
    L.exit        = exit;
    
    ensureTwoVerticalConnectors(L, lotW, lotL);
    finalizeLayout(L, lotW, lotL);
    
    if (hasTooNarrowLane(L)) {
        L.valid = false;
    }
    
    return L;
}


// --------------------------------------------------------
// Horizontal snake generator: just rotate the vertical
// layout by swapping x/y axes and edges.
// --------------------------------------------------------

Edge rotateEdge90(Edge e)
{
    switch (e) {
        case Edge::Bottom: return Edge::Left;
        case Edge::Top:    return Edge::Right;
        case Edge::Left:   return Edge::Bottom;
        case Edge::Right:  return Edge::Top;
    }
    return Edge::Bottom;
}

Layout buildHorizontalSnake(double lotW, double lotL, int N)
{
    // Build as if lot dimensions were swapped (W' = lotL, L' = lotW)
    Layout base = buildVerticalSnake(lotL, lotW, N);
    if (!base.valid) return base;

    Layout L;
    L.valid = true;

    auto rotRect = [](const Rect& r) -> Rect {
        return Rect{ r.y0, r.x0, r.y1, r.x1 };
    };

    L.stalls.reserve(base.stalls.size());
    for (const auto& s : base.stalls)
        L.stalls.push_back({ rotRect(s.r) });

    L.lanes.reserve(base.lanes.size());
    for (const auto& ln : base.lanes)
        L.lanes.push_back({ rotRect(ln.r) });

    L.entrance = base.entrance;
    L.exit     = base.exit;

    L.entrance.edge = rotateEdge90(L.entrance.edge);
    L.exit.edge     = rotateEdge90(L.exit.edge);

    // clean overlaps and compute areas for the real lot orientation
    ensureTwoVerticalConnectors(L, lotW, lotL);
    finalizeLayout(L, lotW, lotL);

    if (hasTooNarrowLane(L)) {
        L.valid = false;
    }

    return L;
}

// --------------------------------------------------------
// Scoring and global search over N + orientation
// --------------------------------------------------------

Layout buildSingleRowLayout(double lotW, double lotL)
{
    Layout L;
    L.valid = false;

    const double lotArea = lotW * lotL;

    // Option A: stalls column along Y, depth along X
    //   - requires: lotW >= STALL_DEPTH and lotL >= STALL_WIDTH
    //   - stalls are 5m (x) by 2.5m (y), repeated along y
    if (lotW >= STALL_DEPTH && lotL >= STALL_WIDTH) {
        const int n = static_cast<int>(std::floor(lotL / STALL_WIDTH));
        if (n > 0) {
            const double x0 = 0.5 * (lotW - STALL_DEPTH);
            const double x1 = x0 + STALL_DEPTH;

            L.stalls.reserve(n);
            for (int i = 0; i < n; ++i) {
                double y0 = i * STALL_WIDTH;
                double y1 = y0 + STALL_WIDTH;
                L.stalls.push_back({ Rect{x0, y0, x1, y1} });
            }

            L.lanes.clear(); // no explicit lane in this special case

            // Simple gates: bottom center -> top center
            L.entrance.edge   = Edge::Bottom;
            L.entrance.width  = STALL_WIDTH;
            L.entrance.offset = 0.5 * (lotW - L.entrance.width);

            L.exit.edge       = Edge::Top;
            L.exit.width      = STALL_WIDTH;
            L.exit.offset     = 0.5 * (lotW - L.exit.width);

            double stallArea = n * STALL_WIDTH * STALL_DEPTH;
            L.stallCount  = n;
            L.usedArea    = stallArea;
            L.wastedArea  = lotArea - stallArea;
            L.valid       = true;
            return L;
        }
    }

    // Option B: rotate 90° – stalls row along X, depth along Y
    //   - requires: lotL >= STALL_DEPTH and lotW >= STALL_WIDTH
    if (lotL >= STALL_DEPTH && lotW >= STALL_WIDTH) {
        const int n = static_cast<int>(std::floor(lotW / STALL_WIDTH));
        if (n > 0) {
            const double y0 = 0.5 * (lotL - STALL_DEPTH);
            const double y1 = y0 + STALL_DEPTH;

            L.stalls.reserve(n);
            for (int i = 0; i < n; ++i) {
                double x0 = i * STALL_WIDTH;
                double x1 = x0 + STALL_WIDTH;
                L.stalls.push_back({ Rect{x0, y0, x1, y1} });
            }

            L.lanes.clear(); // no explicit lane

            // Gates: left center -> right center
            L.entrance.edge   = Edge::Left;
            L.entrance.width  = STALL_WIDTH;
            L.entrance.offset = 0.5 * (lotL - L.entrance.width);

            L.exit.edge       = Edge::Right;
            L.exit.width      = STALL_WIDTH;
            L.exit.offset     = 0.5 * (lotL - L.exit.width);

            double stallArea = n * STALL_WIDTH * STALL_DEPTH;
            L.stallCount  = n;
            L.usedArea    = stallArea;
            L.wastedArea  = lotArea - stallArea;
            L.valid       = true;
            return L;
        }
    }

    // If neither orientation fits even a single stall, this lot is just too small.
    return L;
}

double scoreLayout(const Layout& L)
{
    if (!L.valid) return -1e30;
    // prioritize stall count, then penalize waste
    return L.stallCount * 1000.0 - L.wastedArea;
}

Layout chooseBestLayout(double lotW, double lotL)
{
    Layout best;
    best.valid = false;
    double bestScore = -1e30;

    const double moduleW = 2.0 * STALL_DEPTH + LANE_W;
    int maxN_vert  = static_cast<int>(std::floor(lotW / moduleW));
    int maxN_horiz = static_cast<int>(std::floor(lotL / moduleW));

    // --- vertical snakes ---
    for (int N = 1; N <= maxN_vert; ++N) {
        Layout L = buildVerticalSnake(lotW, lotL, N);
        if (!L.valid) continue;
        double s = scoreLayout(L);
        if (s > bestScore) {
            best      = std::move(L);
            bestScore = s;
        }
    }

    // --- horizontal snakes ---
    for (int N = 1; N <= maxN_horiz; ++N) {
        Layout L = buildHorizontalSnake(lotW, lotL, N);
        if (!L.valid) continue;
        double s = scoreLayout(L);
        if (s > bestScore) {
            best      = std::move(L);
            bestScore = s;
        }
    }

    // --- single-row candidate (even if snakes are valid) ---
    {
        Layout row = buildSingleRowLayout(lotW, lotL);
        if (row.valid) {
            double s = scoreLayout(row);
            if (s > bestScore) {
                best      = std::move(row);
                bestScore = s;
            }
        }
    }

    return best;
}

bool canFitAnyStall(double lotW, double lotL)
{
    // Orientation 1: stall-width along X, stall-depth along Y
    bool orient1 = (lotW >= STALL_WIDTH && lotL >= STALL_DEPTH);

    // Orientation 2: rotated 90° (width ↔ depth)
    bool orient2 = (lotW >= STALL_DEPTH && lotL >= STALL_WIDTH);

    return orient1 || orient2;
}

// --------------------------------------------------------
// SFML drawing
// --------------------------------------------------------

int main()
{
    double lotW, lotL;
    cout << "Enter lot width  (x, metres): ";
    cin  >> lotW;
    cout << "Enter lot length (y, metres): ";
    cin  >> lotL;

    if (!canFitAnyStall(lotW, lotL)) {
        cout << "\nLot is smaller than a single stall in both orientations.\n";
        return 0;
    }

    Layout best = chooseBestLayout(lotW, lotL);

    if (!best.valid) {
        cout << "\nCould not build a snake layout for this lot size.\n";
        return 0;
    }

    cout << "\nBest layout:\n";
    cout << "  Stalls:      " << best.stallCount << "\n";
    cout << "  Used area:   " << best.usedArea   << " m^2\n";
    cout << "  Wasted area: " << best.wastedArea << " m^2\n";
    cout << "  Lanes:       " << best.lanes.size() << "\n";

    // ---------------- SFML window & scaling ----------------

    const float maxWindowWidth  = 1400.0f;
    const float maxWindowHeight = 900.0f;
    const float margin          = 40.0f;

    float scaleX = (maxWindowWidth  - 2.0f * margin) / static_cast<float>(lotW);
    float scaleY = (maxWindowHeight - 2.0f * margin) / static_cast<float>(lotL);
    float scale  = std::min(scaleX, scaleY);
    if (scale <= 0.0f) scale = 1.0f;

    unsigned int windowWidth  =
        static_cast<unsigned int>(lotW * scale + 2.0f * margin);
    unsigned int windowHeight =
        static_cast<unsigned int>(lotL * scale + 2.0f * margin);

    sf::RenderWindow window(
        sf::VideoMode(sf::Vector2u{windowWidth, windowHeight}),
        "Parking Lot Layout"
    );

    auto toScreen = [&](double x, double y) -> sf::Vector2f {
        float sx = margin + static_cast<float>(x) * scale;
        float sy = static_cast<float>(windowHeight) - margin
                 - static_cast<float>(y) * scale;
        return sf::Vector2f(sx, sy);
    };

    sf::RectangleShape lotRect(
        sf::Vector2f(static_cast<float>(lotW) * scale,
                     static_cast<float>(lotL) * scale));
    lotRect.setFillColor(sf::Color(245, 245, 245));
    lotRect.setOutlineThickness(2.0f);
    lotRect.setOutlineColor(sf::Color::Black);
    lotRect.setPosition(toScreen(0.0, lotL));

    vector<sf::RectangleShape> laneShapes;
    laneShapes.reserve(best.lanes.size());
    for (const auto& ln : best.lanes) {
        sf::Vector2f pos = toScreen(ln.r.x0, ln.r.y1);
        sf::RectangleShape r(
            sf::Vector2f(static_cast<float>(ln.r.x1 - ln.r.x0) * scale,
                         static_cast<float>(ln.r.y1 - ln.r.y0) * scale));
        r.setPosition(pos);
        r.setFillColor(sf::Color(190, 190, 190));
        r.setOutlineThickness(1.0f);
        r.setOutlineColor(sf::Color::Black);
        laneShapes.push_back(r);
    }

    vector<sf::RectangleShape> stallShapes;
    stallShapes.reserve(best.stalls.size());
    for (const auto& s : best.stalls) {
        sf::Vector2f pos = toScreen(s.r.x0, s.r.y1);
        sf::RectangleShape r(
            sf::Vector2f(static_cast<float>(s.r.x1 - s.r.x0) * scale,
                         static_cast<float>(s.r.y1 - s.r.y0) * scale));
        r.setPosition(pos);
        r.setFillColor(sf::Color::Transparent);
        r.setOutlineThickness(1.0f);
        r.setOutlineColor(sf::Color(0, 0, 255)); // blueprint-style
        stallShapes.push_back(r);
    }

    while (window.isOpen()) {
        while (const std::optional event = window.pollEvent()) {
            if (event->is<sf::Event::Closed>()) {
                window.close();
            }
        }

        window.clear(sf::Color::White);
        window.draw(lotRect);
        for (const auto& r : laneShapes)   window.draw(r);
        for (const auto& r : stallShapes) window.draw(r);
        window.display();
    }

    return 0;
}
