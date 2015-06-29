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
const double G_TL_SEC = 1e9;
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
            set(x, y, board[y][x] == '.' ? 0 : peg_value[board[y][x] - '0']);
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
            assert(peg(p));
            assert(peg(p + DIFF[dir]));
            assert(!peg(p + 2 * DIFF[dir]));

            sum += at(p + DIFF[dir]);
            move(p, dir);
            p += 2 * DIFF[dir];
        }
        return sum * m.move_dir.size();
    }

    int at(int x, int y) const
    {
        assert(in(x, y));
        return (a[y][x / 16] >> (4 * (x % 16))) & 0xf;
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

    void print() const
    {
        print2d(a, n, n);
    }

    string peg_str() const
    {
        string s;
        rep(y, n)
        {
            rep(x, n)
                s += peg(x, y) ? 'o' : 'x';
            s += '\n';
        }
        return s;
    }

private:
    void set(int x, int y, int peg)
    {
        assert(in(x, y));
        const int shift = 4 * (x % 16);
        a[y][x / 16] = (x % 16 == 15 ? 0 : ((a[y][x / 16] >> (shift + 4)) << (shift + 4))) | (a[y][x / 16] & ((1ull << shift) - 1)) | ((ull)peg << shift);
    }
    void set(const Pos& pos, int peg)
    {
        set(pos.x, pos.y, peg);
    }

    int n;
    ull a[60][4];
};

struct MoveUnit
{
    Pos pos;
    int dir;
    MoveUnit(const Pos& pos, int dir)
        : pos(pos), dir(dir)
    {
    }
    MoveUnit()
#ifndef NDEBUG
        : pos(Pos(-1919810, -114514)), dir(-101019)
#endif
    {
    }
};

template <typename T>
struct Node
{
    T val;
    Node* prev;

    Node()
#ifndef NDEBUG
        :prev(nullptr)
#endif
    {
    }

    Node(const T& val, Node* prev)
        : val(val), prev(prev), size(prev == nullptr ? 1 : prev->size + 1)
    {
    }
    Node(const T& val)
        : val(val), prev(nullptr), size(1)
    {
    }
    Node(Node* prev)
        : prev(prev), size(prev == nullptr ? 1 : prev->size + 1)
    {
    }

    vector<T> list() const
    {
        vector<T> res;
        for (const Node<T>* node = this; node != nullptr; node = node->prev)
            res.push_back(node->val);
        reverse(all(res));
        return res;
    }

    int size;
};
template <typename T>
int size(const Node<T>* node)
{
    return node == nullptr ? 0 : node->size;
}

template <typename T, int SIZE>
class Pool
{
public:
    Pool()
    {
        max_p = 0;
    }

    void init()
    {
        pointer = 0;
    }

    T* get()
    {
        assert(pointer < SIZE);
        upmax(max_p, pointer + 1);
        return &data[pointer++];
    }
    T* get(const T& a)
    {
        assert(pointer < SIZE);
        upmax(max_p, pointer + 1);
        data[pointer] = a;
        return &data[pointer++];
    }


    int max_p;
private:
    T data[SIZE];
    int pointer;
};


class BitBoard
{
public:
    void insert(const Pos& pos)
    {
        f[index(pos)] = true;
    }
    void erase(const Pos& pos)
    {
        f[index(pos)] = false;
    }
    bool count(const Pos& pos) const
    {
        return f[index(pos)];
    }
private:
    int index(const Pos& pos) const
    {
        return (pos.x << 6) | pos.y;
    }
    bitset<64 * 60> f;
};


const int DONE_Q_SIZE = 5;
const int SEARCH_Q_SIZE = 50;

const int BASE_SIZE = (DONE_Q_SIZE + SEARCH_Q_SIZE) * 60 * 60;

Pool<Node<int>, BASE_SIZE * 4> main_move_dir_pool;
Pool<Node<MoveUnit>, BASE_SIZE * 4> move_unit_pool;
Pool<Node<MoveUnit>, BASE_SIZE * 6> move_stack_pool;
Pool<Node<pair<Pos, bool>>, BASE_SIZE * 6 * 4> change_need_stack_pool;

struct State
{
    Board board;
    Pos start_pos;

    BitBoard fixed;
    Pos cur_pos;

    Node<MoveUnit>* move_stack;
    Node<pair<Pos, bool>>* change_need_stack;


    Node<int>* main_move_dir;
    Node<MoveUnit>* prepare_moves;

    State()
        :
            main_move_dir(nullptr),
            prepare_moves(nullptr),
            move_stack(nullptr),
            change_need_stack(nullptr)
    {
    }

    void add_main_move(int dir)
    {
        assert(main_move_dir == nullptr || !board.peg(cur_pos));
        assert(board.peg(cur_pos + DIFF[dir]));
        assert(!board.peg(cur_pos + 2 * DIFF[dir]));

        main_move_dir = main_move_dir_pool.get(Node<int>(dir, main_move_dir));
        fixed.insert(cur_pos + DIFF[dir]);
        fixed.insert(cur_pos + 2 * DIFF[dir]);
        cur_pos += 2 * DIFF[dir]; 
    }

    bool no_prepare_search() const
    {
        return size(move_stack) == 0;
    }

    vector<State> next_main_move_states() const
    {
        assert(no_prepare_search());
        assert(size(move_stack) == 0);
        assert(size(change_need_stack) == 0);
        assert(main_move_dir == nullptr || !board.peg(cur_pos));

        vector<State> next_states;

        rep(dir, 4)
        {
            Pos next = cur_pos + DIFF[dir];
            Pos nextnext = cur_pos + 2 * DIFF[dir];
            if (board.in(nextnext) && !fixed.count(next))
            {
                State nstate = *this;
                nstate.move_stack = move_stack_pool.get(Node<MoveUnit>(MoveUnit(cur_pos, dir), nstate.move_stack));

                // REVIEW: スタック順が両パターンいるか？

                if (!board.peg(nextnext))
                    nstate.fixed.insert(nextnext);
                else
                    nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(nextnext, false), nstate.change_need_stack));

                if (board.peg(next))
                    nstate.fixed.insert(next);
                else
                    nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(next, true), nstate.change_need_stack));

                next_states.push_back(nstate);
            }
        }

        return next_states;
    }

    vector<State> next_states_for_change(const Pos& pos) const
    {
        if (fixed.count(pos))
            return {};

        vector<State> next_states;

        if (board.peg(pos))
        {
            // o -> x

            // oox -> xxo
            // ^ab
            rep(dir, 4)
            {
                Pos a = pos + DIFF[dir];
                Pos b = pos + 2 * DIFF[dir];
                if (board.in(b) && !fixed.count(a) && !fixed.count(b))
                {
                    State nstate = *this;
                    nstate.fixed.insert(pos);
                    nstate.move_stack = move_stack_pool.get(Node<MoveUnit>(MoveUnit(pos, dir), nstate.move_stack));

                    if (board.peg(a))
                        nstate.fixed.insert(a);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(a, true), nstate.change_need_stack));

                    if (!board.peg(b))
                        nstate.fixed.insert(b);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(b, false), nstate.change_need_stack));

                    next_states.push_back(nstate);
                }
            }

            // oox -> xxo
            // a^b
            rep(dir, 4)
            {
                Pos a = pos - DIFF[dir];
                Pos b = pos + DIFF[dir];
                if (board.in(a) && board.in(b) && !fixed.count(a) && !fixed.count(b))
                {
                    State nstate = *this;
                    nstate.fixed.insert(pos);
                    nstate.move_stack = move_stack_pool.get(Node<MoveUnit>(MoveUnit(a, dir), nstate.move_stack));

                    if (board.peg(a))
                        nstate.fixed.insert(a);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(a, true), nstate.change_need_stack));

                    if (!board.peg(b))
                        nstate.fixed.insert(b);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(b, false), nstate.change_need_stack));

                    next_states.push_back(nstate);
                }
            }
        }
        else
        {
            // x -> o

            // xoo -> oxx
            // ^ab
            rep(dir, 4)
            {
                Pos a = pos + DIFF[dir];
                Pos b = pos + 2 * DIFF[dir];
                if (board.in(b) && !fixed.count(a) && !fixed.count(b))
                {
                    State nstate = *this;
                    nstate.fixed.insert(pos);
                    nstate.move_stack = move_stack_pool.get(Node<MoveUnit>(MoveUnit(b, (dir + 2) % 4), nstate.move_stack));

                    if (board.peg(a))
                        nstate.fixed.insert(a);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(a, true), nstate.change_need_stack));

                    if (board.peg(b))
                        nstate.fixed.insert(b);
                    else
                        nstate.change_need_stack = change_need_stack_pool.get(Node<pair<Pos, bool>>(make_pair(b, true), nstate.change_need_stack));

                    next_states.push_back(nstate);
                }
            }
        }

        return next_states;
    }

    vector<State> next_states() const
    {
        if (no_prepare_search())
            return next_main_move_states();
        else
            return next_states_for_change(change_need_stack->val.first);
    }

    void pop_stack()
    {
        while (size(change_need_stack) && board.peg(change_need_stack->val.first) == change_need_stack->val.second)
            change_need_stack = change_need_stack->prev;

        while (size(move_stack) > 1 && board.can_move(move_stack->val.pos, move_stack->val.dir))
        {
            Pos p = move_stack->val.pos;
            int dir = move_stack->val.dir;
            move_stack = move_stack->prev;

            board.move(p, dir);
            assert(board.peg(start_pos));
            prepare_moves = move_unit_pool.get(Node<MoveUnit>(MoveUnit(p, dir), prepare_moves));

            // REVIEW: あやしい
            rep(i, 3)
                fixed.erase(p + i * DIFF[dir]);
            rep(i, 3)
                fixed.insert(move_stack->val.pos + i * DIFF[move_stack->val.dir]);

            while (size(change_need_stack) && board.peg(change_need_stack->val.first) == change_need_stack->val.second)
                change_need_stack = change_need_stack->prev;
        }

        if (size(change_need_stack) == 0 && size(move_stack) == 1)
        {
            Pos p = move_stack->val.pos;
            int dir = move_stack->val.dir;
            move_stack = move_stack->prev;

            assert(main_move_dir == nullptr || !board.peg(p));
            assert(board.peg(p + DIFF[dir]));
            assert(!board.peg(p + 2 * DIFF[dir]));

            add_main_move(dir);
        }
    }

    ll score() const
    {
        return size(main_move_dir);
    }

    bool operator<(const State& other) const
    {
        return score() > other.score();
    }

    vector<Move> make_prepare_moves() const
    {
        vector<Move> moves;
        for (Node<MoveUnit>* node = prepare_moves; node != nullptr; node = node->prev)
            moves.push_back(Move(node->val.pos, {node->val.dir}));
        reverse(all(moves));
        return moves;
    }

    Move make_main_move() const
    {
        return Move(start_pos, main_move_dir->list());
    }
};


vector<Move> search_move(const Board& start_board, const Pos& start)
{
    assert(start_board.at(start));

    move_unit_pool.init();
    main_move_dir_pool.init();
    move_stack_pool.init();
    change_need_stack_pool.init();

    static vector<State> done_q[64 * 64];
    static vector<State> search_q[64 * 64];

    done_q[0].clear();
    search_q[0].clear();

    State start_state;
    start_state.board = start_board;
    start_state.cur_pos = start_state.start_pos = start;
    start_state.fixed.insert(start);
    done_q[0].push_back(start_state);

    int end_qi = 0;
    rep(qi, start_board.size() * start_board.size())
    {
        end_qi = qi;

        if (g_timer.get_elapsed() > G_TL_SEC * 0.9)
            break;

        done_q[qi + 1].clear();
        search_q[qi + 1].clear();

        rep(si, search_q[qi].size())
        {
            State& state = search_q[qi][si];
            state.pop_stack();

            if (state.no_prepare_search())
            {
                done_q[qi].push_back(state);
            }
            else
            {
                if (size(state.move_stack) < 3)
                {
                    vector<State> next_states = state.next_states();

                    auto& q = search_q[qi + 1];
                    for (auto& s : next_states)
                    {
                        if (q.size() < SEARCH_Q_SIZE || s.score() > q.front().score())
                        {
                            q.push_back(s);
                            push_heap(all(q));
                            if (q.size() > SEARCH_Q_SIZE)
                            {
                                pop_heap(all(q));
                                q.pop_back();
                            }
                        }
                    }
                }
            }
        }

        for (const State& state : done_q[qi])
        {
            assert(state.no_prepare_search());
            assert(size(state.move_stack) == 0);

            vector<State> next_states = state.next_states();
            auto& q = search_q[qi + 1];
            for (auto& s : next_states)
            {
                if (q.size() < SEARCH_Q_SIZE || s.score() > q.front().score())
                {
                    q.push_back(s);
                    push_heap(all(q));
                    if (q.size() > SEARCH_Q_SIZE)
                    {
                        pop_heap(all(q));
                        q.pop_back();
                    }
                }
            }
        }

    }

    int best_score = 0;
    vector<Move> best_moves;
    rep(qi, end_qi)
    {
        for (const State& state : done_q[qi])
        {
            assert(state.no_prepare_search());

            vector<Move> prepare_moves = state.make_prepare_moves();

            int score = 0;
            Board board = start_board;
            for (auto& move : prepare_moves)
                score += board.move_score(move);

            Move main_move = state.make_main_move();
            score += board.move_score(main_move);
            if (score > best_score)
            {
                best_score = score;
                best_moves = prepare_moves;
                best_moves.push_back(main_move);
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
        int best_score = 0;
        vector<Move> best;
        rep(y, board.size()) rep(x, board.size())
        {
            if (g_timer.get_elapsed() > G_TL_SEC * 0.95)
                goto TLE;

            if (board.at(x, y))
            {
                auto moves = search_move(board, Pos(x, y));
                if (!moves.empty())
                {
                    Board b = board;
                    int s = 0;
                    for (auto& move : moves)
                        s += b.move_score(move);
                    if (s > best_score)
                    {
                        best_score = s;
                        best = moves;
                    }
                }
            }
        }
TLE:
        if (best.empty())
            break;

        int s = 0;
        for (auto& move : best)
            s += board.move_score(move);
        score += s;
//         fprintf(stderr, "%9d (+%9d)\n", score, s);

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
        dump(main_move_dir_pool.max_p);
        dump(move_unit_pool.max_p);
        dump(move_stack_pool.max_p);
        dump(change_need_stack_pool.max_p);
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
