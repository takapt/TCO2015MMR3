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

// #define USE_RDTSC
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
const double G_TL_SEC = 20;
#else
const double G_TL_SEC = 15.0;
#endif
Timer g_timer;

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

    vector<Pos> make_path() const
    {
        vector<Pos> path;
        Pos p = start;
        path.push_back(p);
        for (int dir : move_dir)
        {
            p += 2 * DIFF[dir];
            path.push_back(p);
        }
        return path;
    }

    Move reverse() const
    {
        Pos goal = start;
        for (int dir : move_dir)
            goal += 2 * DIFF[dir];

        vector<int> rev_dir;
        for (int i = (int)move_dir.size() - 1; i >= 0; --i)
            rev_dir.push_back((move_dir[i] + 2) % 4);

        return Move(goal, rev_dir);
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


int BOARD_SIZE;
class BitBoard
{
public:
    BitBoard()
    {
        clr(f, 0);
    }

    void insert(int x, int y)
    {
        f[y] |= 1ull << x;
    }
    void insert(const Pos& pos)
    {
        insert(pos.x, pos.y);
    }

    void erase(int x, int y)
    {
        f[y] &= ~(1ull << x);
    }
    void erase(const Pos& pos)
    {
        erase(pos.x, pos.y);
    }

    bool count(int x, int y) const
    {
        return f[y] >> x & 1;
    }
    bool count(const Pos& pos) const
    {
        return count(pos.x, pos.y);
    }

    bool peg(int x, int y) const
    {
        return count(x, y);
    }
    bool peg(const Pos& pos) const
    {
        return count(pos);
    }

    bool at(int x, int y) const
    {
        return count(x, y);
    }
    bool at(const Pos& pos) const
    {
        return at(pos);
    }

    bool in(int x, int y) const
    {
        return 0 <= x && x < size() && 0 <= y && y < size();
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

    void move(int x, int y, int dir)
    {
        assert(can_move(x, y, dir));
        insert(x + 2 * DX[dir], y + 2 * DY[dir]);
        erase(x, y);
        erase(x + DX[dir], y + DY[dir]);
    }
    void move(const Pos& pos, int dir)
    {
        move(pos.x, pos.y, dir);
    }

    int size() const
    {
        return ::BOARD_SIZE;
    }

private:
    ull f[60];
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

    int move_score(const MoveUnit& m)
    {
        Pos p = m.pos;
        int dir = m.dir;
        assert(peg(p));
        assert(peg(p + DIFF[dir]));
        assert(!peg(p + 2 * DIFF[dir]));

        int s = at(p + DIFF[dir]);
        move(p, dir);
        p += 2 * DIFF[dir];
        return s;
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

    BitBoard make_bitboard() const
    {
        BitBoard bb;
        rep(y, n) rep(x, n)
            if (peg(x, y))
                bb.insert(x, y);
        return bb;
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
    }

    void init()
    {
        pointer = 0;
    }

    T* get()
    {
        assert(pointer < SIZE);
        return &data[pointer++];
    }
    T* get(const T& a)
    {
        assert(pointer < SIZE);
        data[pointer] = a;
        return &data[pointer++];
    }


private:
    T data[SIZE];
    int pointer;
};




const int SEARCH_DEPTH = 10;
const int DONE_Q_SIZE = 50;
const int SEARCH_Q_SIZE = 150;

const int BASE_SIZE = (DONE_Q_SIZE + SEARCH_Q_SIZE) * 60 * 60;

const int STATE_POOL_SIZE = (DONE_Q_SIZE + SEARCH_Q_SIZE) * 10;

Pool<Node<int>, BASE_SIZE * 4> main_move_dir_pool;
Pool<Node<MoveUnit>, BASE_SIZE * 4> move_unit_pool;
Pool<Node<MoveUnit>, BASE_SIZE * 6> move_stack_pool;
Pool<Node<pair<Pos, bool>>, BASE_SIZE * 6> change_need_stack_pool;

struct State
{
    BitBoard board;
    Pos start_pos;

    BitBoard fixed;
    Pos cur_pos;

    Node<int>* main_move_dir;
    Node<MoveUnit>* prepare_moves;

    Node<MoveUnit>* move_stack;
    Node<pair<Pos, bool>>* change_need_stack;

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


    struct StateFastPrioityQueue
    {
    public:
        StateFastPrioityQueue()
        {
            clr(q_index, -1);
            clear();
        }

        void push(State* state)
        {
            const ll score = state->score();
            assert(score < Q_INDEX_SIZE);
            if (q_index[score] == -1)
            {
                assert(pushed_score.size() < Q_NUM);
                q_index[score] = pushed_score.size();
                pushed_score.push_back(score);
            }
            assert(count(all(pushed_score), score));
            q[q_index[score]].push_back(state);
            ++size_;
        }

        void prepare_pop()
        {
            sort(all(pushed_score));
        }

        State* pop()
        {
            assert(!empty());
            assert(!pushed_score.empty());
            assert(q_index[pushed_score.back()] != -1);

            auto& cur_q = q[q_index[pushed_score.back()]];
            assert(!cur_q.empty());
            State* res = cur_q.back();
            cur_q.pop_back();
            if (cur_q.empty())
            {
                q_index[pushed_score.back()] = -1;
                pushed_score.pop_back();
            }
            --size_;
        }

        void clear()
        {
            for (int i : pushed_score)
            {
                q[q_index[i]].clear();
                q_index[i] = -1;
            }
            pushed_score.clear();
            size_ = 0;
        }

        int size() const
        {
            return size_;
        }
        bool empty() const
        {
            return size() == 0;
        }

    private:
        static const int Q_NUM = 1024; // diffrent scores
        static const int Q_INDEX_SIZE = 2 * ten(6); // need max(score)

        vector<State*> q[Q_NUM];
        int q_index[Q_INDEX_SIZE];
        vector<ll> pushed_score;
        int size_;
    };

    void next_main_move_states(Pool<State, STATE_POOL_SIZE>& state_pool, StateFastPrioityQueue& q) const
    {
        assert(no_prepare_search());
        assert(size(move_stack) == 0);
        assert(size(change_need_stack) == 0);
        assert(main_move_dir == nullptr || !board.peg(cur_pos));

        rep(dir, 4)
        {
            Pos next = cur_pos + DIFF[dir];
            Pos nextnext = cur_pos + 2 * DIFF[dir];
            if (board.in(nextnext) && !fixed.count(next))
            {
                State* p_nstate = state_pool.get(*this);
                State& nstate = *p_nstate;

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

//                 push(q, p_nstate);
                q.push(p_nstate);
            }
        }
    }

    void next_states_for_change(Pool<State, STATE_POOL_SIZE>& state_pool, StateFastPrioityQueue& q) const
    {
        assert(!no_prepare_search());

        const Pos& pos = change_need_stack->val.first;
        if (fixed.count(pos))
            return;

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
                    State* p_nstate = state_pool.get(*this);
                    State& nstate = *p_nstate;

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

//                     push(q, p_nstate);
                    q.push(p_nstate);
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
                    State* p_nstate = state_pool.get(*this);
                    State& nstate = *p_nstate;

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

//                     push(q, p_nstate);
                    q.push(p_nstate);
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
                    State* p_nstate = state_pool.get(*this);
                    State& nstate = *p_nstate;

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

//                     push(q, p_nstate);
                    q.push(p_nstate);
                }
            }
        }
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
#ifndef NDEBUG
            Pos p = move_stack->val.pos;
#endif

            int dir = move_stack->val.dir;
            move_stack = move_stack->prev;

            assert(main_move_dir == nullptr || !board.peg(p));
            assert(board.peg(p + DIFF[dir]));
            assert(!board.peg(p + 2 * DIFF[dir]));

            add_main_move(dir);
        }
    }

    ll ss;
    ll score() const
    {
        return ss;
        ll s = size(main_move_dir) * size(main_move_dir) - size(change_need_stack) + SEARCH_DEPTH * 2;
        assert(s >= 0);
        return s;
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

    struct StatePointerCmp
    {
        bool operator()(const State* a, const State* b)
        {
            return a->score() > b->score();
        }
    };
    void push(vector<State*>& q, State* s) const
    {
        if (q.size() < SEARCH_Q_SIZE || s->score() > q.front()->score())
        {
            q.push_back(s);
            push_heap(all(q), StatePointerCmp());
            if (q.size() > SEARCH_Q_SIZE)
            {
                pop_heap(all(q), StatePointerCmp());
                q.pop_back();
            }
        }
    }

};


bool shuffle_q = false;
struct SearchResult
{
    vector<Move> prepare_moves;
    Move main_move;
    int score;
};

SearchResult search_move(const Board& start_board, const Pos& start, const int search_depth = SEARCH_DEPTH, const BitBoard& start_fixed = BitBoard(), bool extension = false, const map<pair<Pos, Pos>, int>& dist = {})
{
    move_unit_pool.init();
    main_move_dir_pool.init();
    move_stack_pool.init();
    change_need_stack_pool.init();

    static vector<State*> done_q[2];
    static vector<State*> search_q[2];
//     static State::StateFastPrioityQueue done_q[2];
//     static State::StateFastPrioityQueue search_q[2];
    static Pool<State, STATE_POOL_SIZE> state_pool[2];

    static vector<pair<Node<MoveUnit>*, Node<int>*>> done_moves[64 * 64];

//     static State::StateFastPrioityQueue fast_done_q;
    static State::StateFastPrioityQueue fast_search_q;
//     fast_done_q.clear();
    fast_search_q.clear();

    done_q[0].clear();
    search_q[0].clear();
    state_pool[0].init();
    done_moves[0].clear();

    State* start_state = state_pool[0].get(State());
    start_state->board = start_board.make_bitboard();
    start_state->cur_pos = start_state->start_pos = start;
    start_state->fixed = start_fixed;
    start_state->fixed.insert(start);
    done_q[0].push_back(start_state);

    int end_qi = 0;
    rep(qi, start_board.size() * start_board.size())
    {
        const int cur = qi & 1;
        const int next = cur ^ 1;

//         if (done_q[cur].size() + search_q[cur].size() == 0)
//             break;
        if (done_q[cur].empty() && search_q[cur].empty() && fast_search_q.empty())
            break;
        if (qi == 20)
            exit(0);

        if (g_timer.get_elapsed() > G_TL_SEC * 0.95)
            break;

        done_q[next].clear();
        search_q[next].clear();
        state_pool[next].init();
        done_moves[qi].clear();

        if (shuffle_q)
            random_shuffle(all(search_q[cur]));
        for (State* state : search_q[cur])
        {
            state->pop_stack();

            if (state->no_prepare_search())
            {
                done_q[cur].push_back(state);
            }
            else
            {
                if (size(state->move_stack) < search_depth)
                {
//                     state->next_states_for_change(state_pool[next], search_q[next]);
                    state->next_states_for_change(state_pool[next], fast_search_q);
                }
            }
        }

        for (const State* state : done_q[cur])
        {
            end_qi = qi + 1;

            assert(state->no_prepare_search());
            assert(size(state->move_stack) == 0);

            assert(start_fixed.count(state->cur_pos) == dist.count(make_pair(start, state->cur_pos)));
            if (!extension || start_fixed.count(state->cur_pos) && size(state->main_move_dir) > dist.lower_bound(make_pair(start, state->cur_pos))->second)
            {
                done_moves[qi].push_back(make_pair(state->prepare_moves, state->main_move_dir));

//                 if (extension)// && dist[make_pair(start, state->cur_pos)])
//                     fprintf(stderr, "%3d, %3d\n", size(state->main_move_dir), dist[make_pair(start, state->cur_pos)]);
            }

            state->next_main_move_states(state_pool[next], fast_search_q);
        }

        fprintf(stderr, "%4d %3d %3d %3d\n", qi, (int)done_q[cur].size(), (int)search_q[cur].size(), fast_search_q.size());

        fast_search_q.prepare_pop();
        for (int i = 0; i < SEARCH_Q_SIZE && !fast_search_q.empty(); ++i)
            search_q[next].push_back(fast_search_q.pop());
        fast_search_q.clear();
    }

    if (extension)
    {
        SearchResult best;
        best.score = 0;
        rep(qi, end_qi)
        {
            for (auto& it : done_moves[qi])
            {
                static vector<Move> prepare_moves;
                prepare_moves.clear();
                for (auto node = it.first; node != nullptr; node = node->prev)
                    prepare_moves.push_back(Move(node->val.pos, {node->val.dir}));
                reverse(all(prepare_moves));

                static Move main_move;
                main_move.start = start;
                main_move.move_dir.clear();
                for (auto node = it.second; node != nullptr; node = node->prev)
                    main_move.move_dir.push_back(node->val);
                reverse(all(main_move.move_dir));

                Board b = start_board;
                for (auto& move : prepare_moves)
                    b.move(move);

                int sum_peg = 0;
                Pos pos = start;
                for (int dir : main_move.move_dir)
                {
                    sum_peg += b.at(pos + DIFF[dir]);
                    pos += 2 * DIFF[dir];
                }

//                 int score = main_move.move_dir.size() * 10000 + sum_peg;
                int score = ((int)main_move.move_dir.size() - dist.lower_bound(make_pair(start, pos))->second) * 10000 + sum_peg;
                assert(score > 0);
//                 if (pos == start && score > best.score)
                if (score > best.score)
                {
                    assert(!main_move.move_dir.empty());
                    best.score = score;
                    best.prepare_moves = prepare_moves;
                    best.main_move = main_move;
                }
            }
        }
        return best;
    }

    SearchResult best;
    best.score = 0;
    for (int qi = end_qi - 1; qi > max(0, end_qi - 5); --qi)
    {
        for (auto& it : done_moves[qi])
        {
            static vector<Move> prepare_moves;
            prepare_moves.clear();
            for (auto node = it.first; node != nullptr; node = node->prev)
                prepare_moves.push_back(Move(node->val.pos, {node->val.dir}));
            reverse(all(prepare_moves));

            int score = 0;
            Board board = start_board;
            for (auto& move : prepare_moves)
                score += board.move_score(move);

            static Move main_move;
            main_move.start = start;
            main_move.move_dir.clear();
            for (auto node = it.second; node != nullptr; node = node->prev)
                main_move.move_dir.push_back(node->val);
            reverse(all(main_move.move_dir));

            score += board.move_score(main_move);
            if (score > best.score)
            {
                best.score = score;
                best.prepare_moves = prepare_moves;
                best.main_move = main_move;
            }
        }
    }
    return best;
}

SearchResult extend_move(const Board& start_board, const SearchResult& main_result, BitBoard& skip_extend_start)
{
    Board prepared_board = start_board;
    for (auto& move : main_result.prepare_moves)
        prepared_board.move(move);

    BitBoard fixed;
    vector<Pos> ext_pos;
    {
        Pos p = main_result.main_move.start;
        fixed.insert(p);
        ext_pos.push_back(p);
        for (int dir : main_result.main_move.move_dir)
        {
            p += DIFF[dir];
            fixed.insert(p);
            p += DIFF[dir];
            fixed.insert(p);
            ext_pos.push_back(p);
        }
    }

    map<pair<Pos, Pos>, int> dist;
    map<pair<Pos, Pos>, pair<int, int>> index;
    rep(i, ext_pos.size()) rep(j, ext_pos.size())
    {
        auto key = make_pair(ext_pos[i], ext_pos[j]);
        if (dist.count(key))
            upmin(dist[key], abs(j - i));
        else
            dist[key] = abs(j - i);

        if (dist[key] == abs(j - i))
            index[key] = make_pair(i, j);
    }

    SearchResult best_extend;
    best_extend.score = 0;
    for (auto& p : ext_pos)
    {
        if (!skip_extend_start.count(p))
        {
            SearchResult res = search_move(prepared_board, p, 3, fixed, true, dist);
            if (res.score == 0)
                skip_extend_start.insert(p);
            if (res.score > best_extend.score)
                best_extend = res;
        }
    }
    if (best_extend.score == 0)
        return main_result;


    assert(best_extend.main_move.move_dir.size());

    SearchResult extend_res = main_result;
    extend_res.prepare_moves.insert(extend_res.prepare_moves.end(), all(best_extend.prepare_moves));
    {
        Pos start = best_extend.main_move.start;
        Pos goal = best_extend.main_move.make_path().back();
        auto key = make_pair(start, goal);
        assert(dist.count(key));
        if (index[key].first > index[key].second)
        {
            best_extend.main_move = best_extend.main_move.reverse();
            assert(best_extend.main_move.start == goal);
            swap(start, goal);
            swap(key.first, key.second);
        }

        auto& dirs = extend_res.main_move.move_dir;
        dirs.erase(dirs.begin() + index[key].first, dirs.begin() + index[key].second);
        dirs.insert(dirs.begin() + index[key].first, all(best_extend.main_move.move_dir));

//         fprintf(stderr, "%3d, %3d - %3d\n", (int)best_extend.main_move.move_dir.size() - dist[key], (int)best_extend.main_move.move_dir.size(), dist[key]);
    }
//     return main_result;
//     {
//         Pos p = main_result.main_move.start;
//         vector<int> move_dir;
//         bool inserted = false;
//         if (best_extend.main_move.start == p)
//         {
//             inserted = true;
//             move_dir = best_extend.main_move.move_dir;
//         }
//         for (int dir : main_result.main_move.move_dir)
//         {
//             p += 2 * DIFF[dir];
//             move_dir.push_back(dir);
//             if (!inserted && p == best_extend.main_move.start)
//             {
//                 inserted = true;
//                 move_dir.insert(move_dir.end(), all(best_extend.main_move.move_dir));
//             }
//         }
//         assert(inserted);
//         extend_res.main_move.move_dir = move_dir;
//     }

    Board board = start_board;
    int score = 0;
    for (auto& move : extend_res.prepare_moves)
        score += board.move_score(move);
    score += board.move_score(extend_res.main_move);
    extend_res.score = score;

    return extend_res;
}
vector<Move> solve(Board board)
{
    ::BOARD_SIZE = board.size();

    int score = 0;
    vector<Move> res_moves;
    for (int main_move_i = 0; g_timer.get_elapsed() < G_TL_SEC * 0.95; ++main_move_i)
    {
        Timer main_move_timer;
        main_move_timer.start();

        SearchResult best;
        best.score = 0;
        rep(y, board.size()) rep(x, board.size())
//         int x = 2, y = 0;
        {
            const auto skip = [&]()
            {
                if (main_move_i == 0 && g_timer.get_elapsed() > G_TL_SEC * 0.85)
                    return true;
                else if (main_move_i > 0 && g_timer.get_elapsed() > G_TL_SEC * 0.9 && main_move_timer.get_elapsed() > G_TL_SEC * 0.003)
                    return true;
                else if (g_timer.get_elapsed() > G_TL_SEC * 0.95)
                    return true;
                else
                    return false;
            };

            if (skip())
                goto TLE;
//             if (main_move_i == 0 && g_timer.get_elapsed() > G_TL_SEC * 0.9)
//                 goto TLE;
//             else if (main_move_i > 0 && g_timer.get_elapsed() > G_TL_SEC * 0.9 && main_move_timer.get_elapsed() > G_TL_SEC * 0.02)
//                 goto TLE;
//
//             if (g_timer.get_elapsed() > G_TL_SEC * 0.95)
//                 goto TLE;

            if (board.peg(x, y))
            {
                SearchResult res = search_move(board, Pos(x, y));
                BitBoard skip_extend_start;
                while (!skip())
                {
                    SearchResult extend_res = extend_move(board, res, skip_extend_start);
                    if (extend_res.score == res.score)
                        break;
//                     dump((extend_res.main_move.move_dir.size() - res.main_move.move_dir.size()));
//                     fprintf(stderr, "(%2d, %2d): %7d -> %7d\n", x, y, res.score, extend_res.score);
                    res = extend_res;
                }
                if (res.score > best.score)
                {
                    best = res;
//                     fprintf(stderr, "(%2d, %2d): %7d, %4d, %4.1f\n", x, y, best.score, (int)best.main_move.move_dir.size(), g_timer.get_elapsed());
                }
//                     fprintf(stderr, "(%2d, %2d): %7d, %4d, %4.1f\n", x, y, res.score, (int)res.main_move.move_dir.size(), g_timer.get_elapsed());
            }
        }
TLE:
        if (best.score == 0)
            break;

        int s = 0;
        for (auto& move : best.prepare_moves)
            s += board.move_score(move);
        s += board.move_score(best.main_move);
        score += s;
//         fprintf(stderr, "%9d (+%9d)\n", score, s);

        res_moves.insert(res_moves.end(), all(best.prepare_moves));
        res_moves.push_back(best.main_move);
//         break;
    }
    return res_moves;
}

vector<Move> solve_retry(const Board& board)
{
    int best_score = 0;
    vector<Move> best_moves;
    while (g_timer.get_elapsed() < G_TL_SEC * 0.95)
    {
        auto moves = solve(board);
        Board b = board;
        int score = 0;
        for (auto& move : moves)
            score += b.move_score(move);
        if (score > best_score)
        {
            best_score = score;
            best_moves = moves;
        }
        shuffle_q = true;
    }
    return best_moves;
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

void test()
{
    State::StateFastPrioityQueue q;
    State a[5];
    ll s[] = { 3, 1, 4, 1, 5 };
    rep(i, 5)
    {
        a[i].ss = s[i];
        q.push(&a[i]);
    }

    rep(i, 5)
    {
    }
}

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
