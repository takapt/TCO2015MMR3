#ifndef LOCAL
#define NDEBUG
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cerr << a[i]; if (i + 1 != n) cerr << split; } cerr << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cerr.width(width); cerr << a[i][j] << ' '; } cerr << endl; } while (br--) cerr << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)
// #define dump(v)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define erep(i, n) for (int i = 0; i <= (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }

typedef long long ll;
typedef pair<int, int> pint;
typedef unsigned long long ull;

const int DX[] = { 0, 1, 0, -1 };
const int DY[] = { 1, 0, -1, 0 };
const string S_DIR = "DRUL";




ull rdtsc()
{
#ifdef __amd64
    ull a, d;
    __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d)); 
    return (d<<32) | a;
#else
    ull x;
    __asm__ volatile ("rdtsc" : "=A" (x)); 
    return x;
#endif
}
#ifdef LOCAL
const double CYCLES_PER_SEC = 3.30198e9;
#else
const double CYCLES_PER_SEC = 2.5e9;
#endif
double get_absolute_sec()
{
    return (double)rdtsc() / CYCLES_PER_SEC;
}
#ifdef _MSC_VER
#include <Windows.h>
    double get_ms() { return (double)GetTickCount64() / 1000; }
#else
#include <sys/time.h>
    double get_ms() { struct timeval t; gettimeofday(&t, NULL); return (double)t.tv_sec * 1000 + (double)t.tv_usec / 1000; }
#endif

#define USE_RDTSC
class Timer
{
private:
    double start_time;
    double elapsed;

#ifdef USE_RDTSC
    double get_sec() { return get_absolute_sec(); }
#else
    double get_sec() { return get_ms() / 1000; }
#endif

public:
    Timer() {}

    void start() { start_time = get_sec(); }
    double get_elapsed() { return elapsed = get_sec() - start_time; }
};

#ifdef LOCAL
const double G_TL_SEC = 1919810.0;
#else
const double G_TL_SEC = 15.0;
#endif
Timer g_timer;

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int x
             , unsigned int y
             , unsigned int z
             , unsigned int w)
        : x(x), y(y), z(z), w(w) { }
    Random() 
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { assert(upper > 0); return next() % upper; }

    // [low, high]
    int next_int(int low, int high) { assert(low <= high); return next_int(high - low + 1) + low; }

    double next_double(double upper) { return upper * next() / UINT_MAX; }
    double next_double(double low, double high) { return next_double(high - low) + low; }

    template <typename T>
    int select(const vector<T>& ratio)
    {
        T sum = accumulate(ratio.begin(), ratio.end(), (T)0);
        T v = next_double(sum) + (T)1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return 0;
    }
};
Random g_rand;

struct Pos
{
    int x, y;
    Pos(int x, int y)
        : x(x), y(y)
    {
    }
    Pos()
        : x(0), y(0)
    {
    }

    bool operator==(const Pos& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Pos& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Pos& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Pos& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Pos operator+(const Pos& other) const
    {
        Pos res = *this;
        res += other;
        return res;
    }
    Pos operator-(const Pos& other) const
    {
        Pos res = *this;
        res -= other;
        return res;
    }
    Pos operator*(int a) const
    {
        return Pos(x * a, y * a);
    }

    bool operator<(const Pos& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }

    int sq_dist(const Pos& p) const
    {
        int dx = x - p.x;
        int dy = y - p.y;
        return dx * dx + dy * dy;
    }
    double dist(const Pos& p) const
    {
        return sqrt(sq_dist(p));
    }

    bool in_range(const Pos& p, int range) const
    {
        return sq_dist(p) <= range * range;
    }

    Pos next(int dir) const
    {
        return Pos(x + DX[dir], y + DY[dir]);
    }

    void move(int dir)
    {
        x += DX[dir];
        y += DY[dir];
    }
};
Pos operator*(int a, const Pos& pos)
{
    return pos * a;
}
ostream& operator<<(ostream& os, const Pos& pos)
{
    os << "(" << pos.x << ", " << pos.y << ")";
    return os;
}
const Pos DIFF[] = {Pos(DX[0], DY[0]), Pos(DX[1], DY[1]), Pos(DX[2], DY[2]), Pos(DX[3], DY[3])};


struct Move
{
    Pos start;
    vector<int> move_dir;

    Move(){}

    Move(const Pos& start)
        : start(start)
    {
    }

    Move(const Pos& start, const vector<int>& move_dir)
        : start(start), move_dir(move_dir)
    {
    }

    string to_res() const
    {
        stringstream ss;
        assert(!move_dir.empty());
        ss << start.y << ' ' << start.x << ' ';
        for (int dir : move_dir)
            ss << S_DIR[dir];
        return ss.str();
    }
};
class Board
{
public:
    Board(){}
    Board(const vector<int>& peg_value, const vector<string>& board)
        : n(board.size())
    {
        rep(y, n) rep(x, n)
            a[y][x] = board[y][x] == '.' ? 0 : peg_value[board[y][x] - '0'];
    }

    bool in(int x, int y) const
    {
        return 0 <= x && x < n && 0 <= y && y < n;
    }
    bool in(const Pos& pos) const
    {
        return in(pos.x, pos.y);
    }

    bool can_move(int x, int y, int dir) const
    {
        assert(in(x, y));
        return at(x, y)
            && in(x + DX[dir], y + DY[dir]) && at(x + DX[dir], y + DY[dir])
            && in(x + 2 * DX[dir], y + 2 * DY[dir]) && !at(x + 2 * DX[dir], y + 2 * DY[dir]);
    }
    bool can_move(const Pos& pos, int dir) const
    {
        return can_move(pos.x, pos.y, dir);
    }

    bool match(int x, int y, int dir, bool next_peg, bool nextnext_peg) const
    {
        assert(in(x, y));
        return in(x + DX[dir], y + DY[dir]) && in(x + 2 * DX[dir], y + 2 * DY[dir])
            && peg(x + DX[dir], y + DY[dir]) == next_peg && peg(x + 2 * DX[dir], y + 2 * DY[dir]) == nextnext_peg;
    }
    bool match(const Pos& pos, int dir, bool next_peg, bool nextnext_peg) const
    {
        return match(pos.x, pos.y, dir, next_peg, nextnext_peg);
    }
    bool match(int x, int y, int dir, const vector<bool>& peg_pattern) const
    {
        rep(i, peg_pattern.size())
        {
            x += DX[dir];
            y += DY[dir];
            if (!in(x, y) || peg_pattern[i] != peg(x, y))
                return false;
        }
        return true;
    }
    bool match(const Pos& pos, int dir, const vector<bool>& peg_pattern) const
    {
        return match(pos.x, pos.y, dir, peg_pattern);
    }

    void move(int x, int y, int dir)
    {
        assert(can_move(x, y, dir));
        set(x + 2 * DX[dir], y + 2 * DY[dir], at(x, y));
        set(x, y, 0);
        set(x + DX[dir], y + DY[dir], 0);
    }
    void move(const Pos& pos, int dir)
    {
        move(pos.x, pos.y, dir);
    }

    void move(const Move& m)
    {
        Pos p = m.start;
        for (int dir : m.move_dir)
        {
            move(p, dir);
            p += 2 * DIFF[dir];
        }
    }

    int move_score(const Move& m)
    {
        int sum = 0;
        Pos p = m.start;
        for (int dir : m.move_dir)
        {
            sum += at(p + DIFF[dir]);
            move(p, dir);
            p += 2 * DIFF[dir];
        }
        return sum * m.move_dir.size();
    }

    int at(int x, int y) const
    {
        assert(in(x, y));
        return a[y][x];
    }
    int at(const Pos& pos) const
    {
        return at(pos.x, pos.y);
    }

    bool peg(int x, int y) const
    {
        assert(in(x, y));
        return at(x, y);
    }
    bool peg(const Pos& pos) const
    {
        return peg(pos.x, pos.y);
    }

    int size() const
    {
        return n;
    }

private:
    void set(int x, int y, int peg)
    {
        assert(in(x, y));
        a[y][x] = peg;
    }
    void set(const Pos& pos, int peg)
    {
        set(pos.x, pos.y, peg);
    }

    int n;
    int a[64][64];
};


vector<Move> search_move(const Board& start_board, const Pos& start)
{
    assert(start_board.at(start));

    struct State
    {
        Board board;
        Move main_move;
        vector<Move> prepare_moves;

        ull score() const
        {
            ull s = main_move.move_dir.size() * main_move.move_dir.size();
            for (auto& m : prepare_moves)
                s -= m.move_dir.size();
            return s;
        }

        bool operator<(const State& other) const
        {
            return score() > other.score();
        }
    };

    vector<State> q[64 * 64];
    State start_state;
    start_state.board = start_board;
    start_state.main_move.start = start;
    q[0].push_back(start_state);
    rep(qi, start_board.size() * start_board.size())
    {
        sort(all(q[qi]));
        while (q[qi].size() > 20)
            q[qi].pop_back();

        for (const State& state : q[qi])
        {
            set<Pos> cons_peg, cons_empty;
            Pos cur_pos = state.main_move.start;
            cons_peg.insert(cur_pos);
            for (int dir : state.main_move.move_dir)
            {
                cur_pos += DIFF[dir];
                cons_peg.insert(cur_pos);
                cur_pos += DIFF[dir];
                cons_empty.insert(cur_pos);
            }

            rep(main_dir, 4)
            {
                if (state.board.in(cur_pos + 2 * DIFF[main_dir]) && !cons_peg.count(cur_pos + DIFF[main_dir]))
                {
                    if (state.board.match(cur_pos, main_dir, true, false))
                    {
                        State nstate = state;
                        nstate.main_move.move_dir.push_back(main_dir);
                        q[qi + 1].push_back(nstate);
                    }
                    else if (state.board.match(cur_pos, main_dir, {false, true, true})
                            && !cons_peg.count(cur_pos + 2 * DIFF[main_dir])
                            && !cons_peg.count(cur_pos + 3 * DIFF[main_dir])
                            )
                    {
                        State nstate = state;
                        nstate.board.move(cur_pos + 3 * DIFF[main_dir], (main_dir + 2) % 4);
                        nstate.prepare_moves.push_back(Move(cur_pos + 3 * DIFF[main_dir], {(main_dir + 2) % 4}));
                        nstate.main_move.move_dir.push_back(main_dir);
                        q[qi + 2].push_back(nstate);
                    }
                }
            }
        }
    }

    int best_score = 0;
    vector<Move> best_moves;
    rep(qi, start_board.size() * start_board.size())
    {
        for (auto& state : q[qi])
        {
            int score = 0;
            Board board = start_board;
            for (auto& move : state.prepare_moves)
                score += board.move_score(move);
            score += board.move_score(state.main_move);
            if (score > best_score)
            {
                best_score = score;
                best_moves = state.prepare_moves;
                best_moves.push_back(state.main_move);
            }
        }
    }
    return best_moves;
}

vector<Move> solve(Board board)
{
    int score = 0;
    vector<Move> res_moves;
    while (g_timer.get_elapsed() < G_TL_SEC * 0.95)
    {
        vector<Move> best;
        rep(y, board.size()) rep(x, board.size())
        {
            if (g_timer.get_elapsed() > G_TL_SEC * 0.95)
                goto TLE;

            if (board.at(x, y))
            {
                auto moves = search_move(board, Pos(x, y));
                if (!moves.empty() && (best.empty() || moves.back().move_dir.size() > best.back().move_dir.size()))
                    best = moves;
            }
        }
TLE:
        if (best.empty())
            break;

        int s = 0;
        for (auto& move : best)
            s += board.move_score(move);
        score += s;
        fprintf(stderr, "%6d (+%6d)\n", score, s);

        res_moves.insert(res_moves.end(), all(best));
    }
    return res_moves;
}

class PegJumping
{
public:
    vector<string> getMoves(const vector<int>& peg_value, const vector<string>& board)
    {
        g_timer.start();

        vector<string> res;
        for (auto& move : solve(Board(peg_value, board)))
            res.push_back(move.to_res());

        dump(g_timer.get_elapsed());
        return res;
    }
};


#ifdef LOCAL
int main()
{
    int m;
    cin >> m;
    vector<int> pegValue(m);
    input(pegValue, m);

    int n;
    cin >> n;
    vector<string> board(n);
    input(board, n);

    vector<string> ret = PegJumping().getMoves(pegValue, board);
    cout << ret.size() << endl;
    for (auto& s : ret)
        cout << s << endl;
    cout.flush();
}
#endif
