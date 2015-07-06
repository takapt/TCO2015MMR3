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
            if (score >= Q_INDEX_SIZE)
                dump(score);
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

            return res;
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

