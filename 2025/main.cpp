/*
@project : 2025嵌入式大赛
@data : 2025.5.23
*/

#pragma GCC optimize("Ofast,inline,unroll-loops")
#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <climits>
#include <numeric>
#include <cassert>
#include <chrono>
#include <random>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <cstring>
using namespace std;

#define INLINE   inline __attribute__ ((always_inline))
#define NOINLINE __attribute__ ((noinline))
#define h(x) pow(2, -1.0 * x / 100)
#define p(x) pow(2, -1.0 * x / 200)
#define x first
#define y second

#ifndef LOCAL
#define LOCAL
#endif
#ifndef TIMEOUTINFO
// #define TIMEOUTINFO
#endif
#ifndef NPUINFO
// #define NPUINFO
#endif
typedef pair<int, int> PII;

constexpr int LENGTH = 6e5;
constexpr int SEGMENT = 6;
constexpr int MAXM = 500 + 10;
constexpr int  MAXT = 300 + 10;
constexpr int MAXN = 10 + 10;
constexpr int MAXB = 1000 + 10;
constexpr int LIMIT = 0;

/*************************************************************************************************/
const auto startTime = std::chrono::steady_clock::now();

 int runtime() {
    auto now = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - startTime);
    return int(duration.count());
}
/*************************************************************************************************/
/**随机数生成器*/
struct RNG {
	unsigned int x = 123456789;
	unsigned int y = 362436069;
	unsigned int z = 521288629;
    unsigned int rand() {
		x ^= x << 16;
		x ^= x >> 5;
		x ^= x << 1;
		unsigned int t = x;
		x = y; y = z; z = t ^ x ^ y;
		return z;		
    }
    //  int next(int x) {return ((long long)rand() * x) >> 32;}
     int next(int x) {return rand() % x;}
     int next(int a, int b) {return a + (rand() % (b - a));}
    //  int next(int a, int b) {return a + ((long long) rand() * (b - a)) >> 32;}
     double next_double() {return (rand() + 0.5) * (1.0 / 4294967296.0);}
}; 
static RNG rng;

struct Request {
    int arrive_time;  
    double process_time;
    int memory;          
    int user_id;         
    int batch;
    int index;       
    Request(){

    }
    Request(int _arrive_time, int _process_time, int _memory, int _user_id, int _batch, int _index) : 
    arrive_time(_arrive_time), process_time(_process_time), memory(_memory), user_id(_user_id), batch(_batch), index(_index){}

    bool operator<(const Request& other) const {
        if (arrive_time == other.arrive_time)
            return user_id > other.user_id; 
        return arrive_time > other.arrive_time;
    }
};

struct Sample {
    int arriveTime;  
    int processTime;
    int endTime;
    int userId;      
    int batch;     
    int index;       

    Sample(){

    }
    Sample(int _arrive_time, int _process_time, int _end_time, int _user_id, int _batch, int _index) : 
    arriveTime(_arrive_time), processTime(_process_time), endTime(_end_time), userId(_user_id), batch(_batch), index(_index){}

    bool operator<(const Sample& other) const {
        if (arriveTime == other.arriveTime)
            return userId < other.userId; 
        return arriveTime < other.arriveTime;
    }
};
struct NPUState {
    priority_queue<Request> pending; // 等待处理的请求队列
    vector<Sample> waitSamples; // 保存 0: 到达时间; 1 : 用户编号; 2 : batch下标, 3:结束时间
    vector<PII> schedule; // 按时间排序的时间段
    int useMemory[LENGTH];
    // 统计放入后区间段[用户，样本]
    vector<tuple<int, int, int>> infos;               
    int current_time = 0.0;        


    void reset(bool flag) {
        waitSamples.clear();
        current_time = 0;
        infos.clear();
        memset(useMemory, 0, sizeof useMemory);
        if (flag) {
            // 清空优先队列
            while(!pending.empty()) pending.pop();
        }
    }     
};

struct Server {
    int g;       
    int k;       
    int m;       
    vector<NPUState> npus; 
};

struct UserSolution {
    int id{};
    int T{};
    vector<int> sendTimes;
    vector<int> arriveTimes;
    vector<int> endTimes;
    vector<int> servers;
    vector<int> npus;
    vector<int> batches;
    int moveCount = 0;
    int endTime = 0;
    bool valid = false;
};

// 预处理数据部分
double hxInt[MAXM], pxInt[MAXT];
int expense[MAXN][MAXB];

void preProcess() {
    for (int i = 0; i < MAXM; ++ i) {
        hxInt[i] = h(i);
    }
    for (int i = 0; i < MAXT; ++ i) {
        pxInt[i] = p(i);
    }
}

class Scheduler {
private:
    int N, M, a, b;
    vector<Server> servers;
    vector<tuple<int, int, int>> users; 
    vector<vector<int>> latency;        
    vector<UserSolution> bestSolutions;
    vector<bool> TimeOut;
    double bestScore{};

    vector<PII> latencyBound;
    vector<PII> batchBound;

public:
    void init() {
        scanf("%d", &N);
        servers.resize(N);
        for (int i = 0; i < N; ++i) {
            scanf("%d%d%d", &servers[i].g, &servers[i].k, &servers[i].m);
            servers[i].npus.resize(servers[i].g);
        }
        
        scanf("%d", &M);
        users.resize(M);
        for (int i = 0; i < M; ++i) {
            int s, e, cnt;
            scanf("%d%d%d", &s, &e, &cnt);
            users[i] = {s, e, cnt};
        }
        
        latency.resize(N);
        latencyBound.resize(M);
        batchBound.resize(M);
        TimeOut.assign(M, 0);
        for (int i = 0; i < M; ++ i) {
            latencyBound[i] = make_pair(INT_MAX, INT_MIN);
            batchBound[i] = make_pair(INT_MAX, INT_MIN);
        }
        for (int i = 0; i < N; ++i) {
            latency[i].resize(M);
            for (int j = 0; j < M; ++j) {
                scanf("%d", &latency[i][j]);
                latencyBound[j].x = min(latencyBound[j].x, latency[i][j]);
                latencyBound[j].y = max(latencyBound[j].y, latency[i][j]);
            }
        }
        scanf("%d%d",&a, &b);



        // 预处理计算处理时间
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < MAXB; ++j) {
                expense[i][j] = ceil(sqrt(j) / servers[i].k);
            }
        }

        // 预处理用户是否一定会超时
        for (int user_id = 0; user_id < M; ++ user_id ) {
            auto& [s_i, e_i, cnt_i] = users[user_id];
            // 计算批次和到达时间
            int batch = batchBound[user_id].y;
            int T = (cnt_i + batch - 1) / batch;
            double last_arrive = s_i + (T-1)*(latencyBound[user_id].x+1) + latencyBound[user_id].x;
            if (last_arrive >= e_i) {
                TimeOut[user_id] = true;
            }
        }
    }

    void schedulingOne(vector<vector<int>>& minSegment, int segment) {
        // 随机种子

        // 构造初始服务器
        vector<PII> allServerNpu;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < servers[i].g; ++j) {
                allServerNpu.emplace_back(i, j);
            }
        }
        // 生成用户列表
        vector<PII> userGroups(M);
        for (int i = 0; i < M; ++i) {
            userGroups[i] = make_pair(i, get<2>(users[i]));
        }
        sort(userGroups.begin(), userGroups.end(), [&](const PII &p1, const PII &p2) {
            return p1.y < p2.y;
        });
        int timeOutCnt = 0;
        for (int is = 0; is < 3; ++is) {
            if (is == 1) {
                for (int i = 0; i < M; ++i) {
                    userGroups[i] = make_pair(i, get<1>(users[i]));
                }
                sort(userGroups.begin(), userGroups.end(), [&](const PII &p1, const PII &p2) {
                    return p1.y < p2.y;
                });
            } else if (is == 2) {
                for (int i = 0; i < M; ++i) {
                    userGroups[i] = make_pair(i, get<1>(users[i]) *get<2>(users[i]) );
                }
                sort(userGroups.begin(), userGroups.end(), [&](const PII &p1, const PII &p2) {
                    return p1.y < p2.y;
                });      
            }
            int op = 0;
            // 更新服务器信息
            retSetAllNpuState(true);
            vector<UserSolution> solutions(M);
            // 第一阶段：为每个用户生成请求
            int sz = allServerNpu.size();

            vector<int> occupyTime(sz, 0);
            vector<PII> remainUserGroups;
            for (auto [userId, _]: userGroups) {
                auto [s_i, e_i, cnt_i] = users[userId];
                UserSolution sol;
                sol.id = userId;
                int earliestSendTime = s_i;
                int endTime = 0;
                sol.T = 0;
                // 遍历所有npu
                int minOccupyTime = INT_MAX;
                for (int i = 0; i < sz; ++i) {
                    if (occupyTime[i] < minOccupyTime) {
                        minOccupyTime = occupyTime[i];
                        op = i;
                    }
                }
                int bestServer = allServerNpu[op].x, bestNpu = allServerNpu[op].y;
                auto &npu = servers[bestServer].npus[bestNpu];
                vector<int> saveServers, saveNpus, saveBatches, saveSendTimes, saveArriveTime;
                int maxBatch = min((servers[bestServer].m - b) / a, 1000);
                // 调整
                segment = min(segment, minSegment[bestServer][userId]);
                maxBatch = max(1, maxBatch / segment);
                if (maxBatch  * 300 < cnt_i) {
                    return;
                }
                while (cnt_i) {
                    // 定位发送时刻
                    int tempArriveTime = earliestSendTime + latency[bestServer][userId];
                    // 寻找调度区间的起始点
                    int finalBatch = -1, realBatch = -1, finalTime = -1;
                    while (1) {
                        // 此时的tempArriveTime时刻定位在调度区间的起始点（获得此时剩余内存对应的batch)
                        while ((servers[bestServer].m - npu.useMemory[tempArriveTime] - b) / a < 1) {
                            ++tempArriveTime;
                        }
                        if (tempArriveTime > LENGTH / 2) {
                            return;
                        }
                        int remainMemory = servers[bestServer].m - npu.useMemory[tempArriveTime];
                        // 可以调度的batchA
                        int batchA = min((remainMemory - b) / a, 1000);
                        batchA = min(batchA, maxBatch);
                        // 实际调度的batch
                        int realBatchA = min(batchA, cnt_i);
                        int remain = 300 - sol.T - 1;
                        vector<int> canPutBatch(20, 0); // 放置20秒的目前可以调度的batch
                        for (int j = 0; j < 20; ++j) {
                            int curTime = tempArriveTime + j;
                            int remainMemory = servers[bestServer].m - npu.useMemory[curTime];
                            canPutBatch[j] = max(0, min((remainMemory - b) / a, 1000));
                        }
                        int finalBatchA = -1, finalTimeA = -1;
                        // 计算满足的最终最大可以调度的batch(用二分)
                        for (int batch = realBatchA; batch > 0; --batch) {
                            int costTime = expense[bestServer][batch];
                            bool flag = true;
                            for (int j = 0; j < costTime; ++j) {
                                if (canPutBatch[j] < batch) {
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag) {
                                finalBatchA = batch;
                                finalTimeA = (tempArriveTime + costTime - 1);
                                break;
                            }
                        }
                        // 满足条件
                        if (finalBatchA >= 1 && realBatchA * remain + finalBatchA >= cnt_i) {
                            finalBatch = finalBatchA;
                            finalTime = finalTimeA;
                            realBatch = realBatchA;
                            break;
                        }
                        ++tempArriveTime;
                    }
                    // 可以调度
                    finalBatch = min(finalBatch, realBatch);
                    int finalMemory = finalBatch * a + b;
                    int arriveTime = tempArriveTime;
                    int bestSendTime = arriveTime - latency[bestServer][userId];
                    earliestSendTime = bestSendTime;
                    for (int j = tempArriveTime; j <= finalTime; ++j) {
                        if (!npu.useMemory[j]) ++occupyTime[op];
                        npu.useMemory[j] += finalMemory;
                    }
                    saveSendTimes.emplace_back(earliestSendTime);
                    saveArriveTime.emplace_back(arriveTime);
                    saveServers.emplace_back(bestServer + 1);
                    saveNpus.emplace_back(bestNpu + 1);
                    saveBatches.emplace_back(finalBatch);
                    earliestSendTime = bestSendTime + latency[bestServer][userId] + 1;
                    endTime = max(endTime, finalTime);
                    cnt_i -= finalBatch;
                    ++sol.T;
                }
                if (endTime + 1 > e_i) {
                    //撤回
                    for (int i = 0; i < saveServers.size(); ++i) {
                        int accupy = expense[saveServers[i] - 1][saveBatches[i]];
                        int sub = 150;
                        for (int j = max(1,saveSendTimes[i] - sub); j < max(1,saveSendTimes[i] - sub) + accupy; ++j) {
                            npu.useMemory[j] -= saveBatches[i] * a + b;
                            if (!npu.useMemory[j]) --occupyTime[op];
                        }
                    }
                    remainUserGroups.emplace_back(userId, 0);
                } else {
                    sol.sendTimes = move(saveSendTimes);
                    sol.arriveTimes = move(saveArriveTime);
                    sol.servers = move(saveServers);
                    sol.npus = move(saveNpus);
                    sol.batches = move(saveBatches);
                    sol.endTime = endTime;
                    solutions[userId] = move(sol);
                }
            }
            for (auto [userId, _]: remainUserGroups) {
                auto [s_i, e_i, cnt_i] = users[userId];
                UserSolution sol;
                sol.id = userId;
                int earliestSendTime = s_i;
                int endTime = 0;
                sol.T = 0;
                // 遍历所有npu
                int minOccupyTime = INT_MAX;
                for (int i = 0; i < sz; ++i) {
                    if (occupyTime[i] < minOccupyTime) {
                        minOccupyTime = occupyTime[i];
                        op = i;
                    }
                }
                int bestServer = allServerNpu[op].x, bestNpu = allServerNpu[op].y;
                auto &npu = servers[bestServer].npus[bestNpu];
                vector<int> saveServers, saveNpus, saveBatches, saveSendTimes, saveArriveTime;
                int maxBatch = min((servers[bestServer].m - b) / a, 1000);
                // 调整
                maxBatch = max(1, maxBatch / segment);
                if (maxBatch * 300 < cnt_i) {
                    return;
                }
                while (cnt_i) {
                    // 定位发送时刻
                    int tempArriveTime = earliestSendTime + latency[bestServer][userId];
                    // 寻找调度区间的起始点
                    int finalBatch = -1, realBatch = -1, finalTime = -1;
                    while (1) {
                        // 此时的tempArriveTime时刻定位在调度区间的起始点（获得此时剩余内存对应的batch)
                        while ((servers[bestServer].m - npu.useMemory[tempArriveTime] - b) / a < 1) {
                            ++tempArriveTime;
                        }
                        if (tempArriveTime > LENGTH / 2) {
                            return;
                        }
                        int remainMemory = servers[bestServer].m - npu.useMemory[tempArriveTime];
                        // 可以调度的batchA
                        int batchA = min((remainMemory - b) / a, 1000);
                        batchA = min(batchA, maxBatch);
                        // 实际调度的batch
                        int realBatchA = min(batchA, cnt_i);
                        int remain = 300 - sol.T - 1;
                        vector<int> canPutBatch(20, 0); // 放置20秒的目前可以调度的batch
                        for (int j = 0; j < 20; ++j) {
                            int curTime = tempArriveTime + j;
                            int remainMemory = servers[bestServer].m - npu.useMemory[curTime];
                            canPutBatch[j] = max(0, min((remainMemory - b) / a, 1000));
                        }
                        int finalBatchA = -1, finalTimeA = -1;
                        // 计算满足的最终最大可以调度的batch(用二分)
                        for (int batch = realBatchA; batch > 0; --batch) {
                            int costTime = expense[bestServer][batch];
                            bool flag = true;
                            for (int j = 0; j < costTime; ++j) {
                                if (canPutBatch[j] < batch) {
                                    flag = false;
                                    break;
                                }
                            }
                            if (flag) {
                                finalBatchA = batch;
                                finalTimeA = (tempArriveTime + costTime - 1);
                                break;
                            }
                        }
                        // 计算
                        // finalBatchA = min(finalBatchA, realBatchA);
                        // 满足条件
                        if (finalBatchA >= 1 && realBatchA * remain + finalBatchA >= cnt_i) {
                            finalBatch = finalBatchA;
                            finalTime = finalTimeA;
                            realBatch = realBatchA;
                            break;
                        }
                        ++tempArriveTime;
                    }
                    // 可以调度
                    finalBatch = min(finalBatch, realBatch);
                    int finalMemory = finalBatch * a + b;
                    int arriveTime = tempArriveTime;
                    int bestSendTime = arriveTime - latency[bestServer][userId];
                    earliestSendTime = bestSendTime;
                    for (int j = tempArriveTime; j <= finalTime; ++j) {
                        if (!npu.useMemory[j]) ++occupyTime[op];
                        npu.useMemory[j] += finalMemory;
                    }
                    saveSendTimes.emplace_back(earliestSendTime);
                    saveArriveTime.emplace_back(arriveTime);
                    saveServers.emplace_back(bestServer + 1);
                    saveNpus.emplace_back(bestNpu + 1);
                    saveBatches.emplace_back(finalBatch);
                    earliestSendTime = bestSendTime + latency[bestServer][userId] + 1;
                    endTime = max(endTime, finalTime);
                    cnt_i -= finalBatch;
                    ++sol.T;
                }

                sol.sendTimes = move(saveSendTimes);
                sol.arriveTimes = move(saveArriveTime);
                sol.servers = move(saveServers);
                sol.npus = move(saveNpus);
                sol.batches = move(saveBatches);
                sol.endTime = endTime;
                solutions[userId] = move(sol);

            }

            function<void(vector<UserSolution>&)> constructQueue = [&](vector<UserSolution>& solutions)->void {
                // 根据解决方案构造队列
                retSetAllNpuState(true);
                for (auto& user : solutions) user.endTime = 0, user.valid = false, user.endTimes.assign(user.batches.size(), 0);
                // 放进vector中{userId, batchIndex}, 所有的npu
                vector<vector<vector<PII>>> allNpuQueue(N);
                for (int s = 0; s < N; ++ s) {
                    allNpuQueue[s].resize(servers[s].g);
                }
                // 将用户的所有batch放入队列中
                for (auto& user : solutions) {
                    int sz = user.batches.size();
                    assert(user.batches.size() == user.T && user.servers.size() == user.T && user.sendTimes.size() == user.T && user.npus.size() == user.T);
                    for (int i = 0; i < sz; ++ i) {
                        int server = user.servers[i] - 1;
                        int npu = user.npus[i] - 1;
                        allNpuQueue[server][npu].emplace_back(user.id, i);
  
                    }
                }

                // 排序所有队列
                for (int s = 0; s < N; ++ s) {
                    for (int npu = 0; npu < servers[s].g; ++ npu) {
                        sort(allNpuQueue[s][npu].begin(), allNpuQueue[s][npu].end(), [&](const PII& task1, const PII& task2){
                            auto& time1 = solutions[task1.x].arriveTimes[task1.y];
                            auto& time2 = solutions[task2.x].arriveTimes[task2.y];
                            if (time1 != time2) {
                                return time1 < time2;
                            }
                            return task1.x < task2.x;
                        });
                    }
                }

                // 对所有用户记忆化存储
                vector<vector<bool>> st(M);
                for (auto& user : solutions) {
                    int sz = user.batches.size();
                    st[user.id].assign(sz, 0);
                }

                // 开始调度
                for (int s = 0; s < N; ++ s) {
                    for (int npu = 0; npu < servers[s].g; ++ npu) {
                        auto& q = allNpuQueue[s][npu];
                        auto& npuMachine = servers[s].npus[npu];
                        int total = q.size();
                        // 遍历
                        int start = 0;
                        while (total) {
                            // 遍历队列
                            for (const auto&[userId, batchIndex] : q) {
                                // 已经调度
                                if (st[userId][batchIndex]) continue;
                                // 调度请求已经大于当前时间
                                if (solutions[userId].arriveTimes[batchIndex] > start) {
                                    break;
                                }
                                // 此时请求满足可以调度，需要判断内存是否满足
                                int needBatch = solutions[userId].batches[batchIndex];
                                int needMemory = a * needBatch + b;
                                // 遍历调度时间
                                int processTime = expense[s][needBatch];
                                int endTime = start + processTime - 1;
                                bool isUse = true;
                                // cerr << start << " " << endTime << endl;
                                for (int i = start; i <= endTime; ++ i) {
                                    if (npuMachine.useMemory[i] + needMemory > servers[s].m) {
                                        isUse = false;
                                        break;
                                    }
                                }
                                if (isUse) {
                                    // 可以调度
                                    for (int i = start; i <= endTime; ++ i) {
                                        npuMachine.useMemory[i] += needMemory;
                                    }
                                    solutions[userId].endTimes[batchIndex] = endTime;
                                    solutions[userId].endTime = max(solutions[userId].endTime, endTime);
                                    st[userId][batchIndex] = true;
                                    -- total;
                                }
                            }
                            ++ start;
                        }

                    }
                }

            };
            // constructQueue(solutions);
            double curScore = getScore(solutions, timeOutCnt);
            cerr << "分段为: " << segment << "下的情况: " << is << endl;
            cerr << "当前分数: " << curScore << endl;
            cerr << "超时的数量: " << timeOutCnt << endl;
            if (curScore > bestScore) {
                bestScore = curScore;
                bestSolutions = solutions;
                cerr << "当前贪心解得分:" << fixed <<  bestScore << endl;
                cerr << "理论最优分数: " << fixed << calculateUpperScore() << endl;
                cerr << "理论最优分数（占比): " << 1.0 * bestScore / calculateUpperScore() << endl;
                cerr << "超时用户: " << timeOutCnt << "; 总用户: " << M << endl;
            }
        }

        // cerr << "best Score: " << getScore(bestSolutions, timeOutCnt) << endl;
        // cerr << "超时的数量: " << timeOutCnt << endl;
    }

    INLINE bool check(const vector<UserSolution> &userSolution) const {

        // 检查每个用户的样本推理时间是不是严格递增并且大于s_i
        for (const auto& sol : userSolution) {
            auto [s_i, e_i, cnt_i] = users[sol.id];
            int cnt = 0;

            if (sol.T < 1 || sol.T > 300) {
                return false;
            }
            for (int i = 0; i < sol.T; ++ i) {
                int serverId = sol.servers[i] - 1, npuId = sol.npus[i] - 1;
                int batch = sol.batches[i];
                int sendTime = sol.sendTimes[i];
                int maxBatch = min((servers[serverId].m - b)/a, 1000);
                if (serverId < 0 || serverId >= N || npuId < 0 || npuId >= servers[serverId].g) return false;
                if (batch > maxBatch || sendTime < s_i +  i*(latency[serverId][sol.id]+1)) return false;
                
                cnt += batch;
            }
            if (cnt_i != cnt) return false;
        }
        return true;
    }
    void schedulingTwo(int startTime) {
        vector<UserSolution> secondBestSolutions(M);
        double secondBestScore = 0;
        // 随机种子
        random_device rd;
        mt19937 gen(1);
        static int iter = 0;
        // 构造初始服务器
        vector<PII> allServerNpu;
        for (int i = 0; i < N; ++ i) {
            for (int j = 0; j < servers[i].g; ++ j) {
                allServerNpu.emplace_back(i, j);
            }
        }
        int selectOp = allServerNpu.size();
        // 生成用户列表
        vector<PII> userGroups(M);
        for (int i = 0; i < M; ++ i) {
            userGroups[i] = make_pair(i, get<2>(users[i]));
        }
        sort(userGroups.begin(), userGroups.end(), [&](const PII& p1, const PII& p2){
            return p1.y < p2.y ;
        });
        for (int is = 0; is < 2; ++ is) {
            
            if (is == 1) {
                for (int i = 0; i < M; ++ i) {
                    userGroups[i] = make_pair(i, get<1>(users[i]));
                }
                sort(userGroups.begin(), userGroups.end(), [&](const PII& p1, const PII& p2){
                    return p1.y < p2.y ;
                });                
            }


            int op = 0;
            // 更新服务器信息
            retSetAllNpuState(true);
            int timeOutCnt = 0;
            vector<UserSolution> solutions(M);

            vector<int> timeOutCandidates;
            // 第一阶段：为每个用户生成请求
            int up = allServerNpu.size();
            for (auto [user_id, _] : userGroups) {
                int num = 0;
                while (num ++ < up) {
                    auto& [s_i, e_i, cnt_i] = users[user_id];
                    UserSolution sol;
                    sol.id = user_id;

                    // 寻找最优服务器和NPU(轮询)
                    int best_server = allServerNpu[op].x, best_npu = allServerNpu[op].y;
                    op = (op + 1) % selectOp;
                    

                    // 生成批次和请求
                    int max_batch = min((servers[best_server].m - b) / a, 1000);
                    int batch = min(max_batch, cnt_i);
                    vector<int> batches;
                    int minTime = findOptimalSchedule(cnt_i, servers[best_server].k, batch, batches);
                    int T = batches.size();

                    int startSendTime = s_i;
                    int fixedLatency = latency[best_server][user_id]+1;
                    // 为每一个批次确定发送时间
                    auto& npu = servers[best_server].npus[best_npu];
                    // 保存npu的旧状态
                    vector<PII> backupSchedule = npu.schedule;
                    vector<Sample> backupWaitSample = npu.waitSamples;
                    vector<int> saveServers(T), saveNpus(T), saveBatches(T), saveSendTimes(T), saveArriveTime(T);
                    bool isTimeOut = false;
                    for (int i = 0; i < T; ++ i) {

                        if (isTimeOut) break;
                        // 计算到达时间：发送时间 + 延迟
                        int earliestSendTime = startSendTime;
                        int arriveTime = earliestSendTime + latency[best_server][user_id];
                        int processTime = expense[best_server][batches[i]];


                        // 在npu等待样本中寻找位置
                        Sample sample{arriveTime, processTime, INT_MAX, user_id, batches[i], i};
                        auto it = upper_bound(npu.waitSamples.begin(), npu.waitSamples.end(), sample);
                        int pos = it - npu.waitSamples.begin();
                        // pos : 表示该位置大于我寻找的位置，所以应该是前一个
                        pos -= 1;
                        int preUserId = -1, preUserBatchId = -1, preEndTime = -1;
                        int canExcTime = -1;
                        if (pos < 0) {
                            // 表示我的插入位置在第一个
                            preEndTime = 1;
                            // 寻找包含1所在区间的末端点
                            for (auto& [l, r] : npu.schedule) {
                                if (l <= preEndTime) {
                                    preEndTime = r;
                                    break;
                                }
                            }
                            canExcTime = preEndTime;


                        } else {
                            preUserId = npu.waitSamples[pos].userId, preUserBatchId = npu.waitSamples[pos].index, preEndTime = npu.waitSamples[pos].endTime;
                            canExcTime = preEndTime;
                            // 找到包含该时间的区间
                            for (auto& [l, r] : npu.schedule) {
                                if (canExcTime >= l && canExcTime <= r) {
                                    canExcTime = r;
                                    break;
                                }
                            }
                        }
                        // // 找到了连续执行区间的末端点, 调整到达时间和发送时间
                        arriveTime = max(arriveTime, canExcTime);
                        earliestSendTime = max(earliestSendTime, arriveTime - latency[best_server][user_id]);
                        startSendTime = earliestSendTime + fixedLatency;

                        // 将[user, batch] 插入到npu.waitSamples中
                        sample.arriveTime = arriveTime;
                        auto newIt = upper_bound(npu.waitSamples.begin(), npu.waitSamples.end(), sample);
                        npu.waitSamples.insert(newIt, sample);
                        int newPos = pos + 1;


                        // 更新waitSamples和schedules
                        npu.current_time = 0;
                        npu.schedule.clear();

                        int l = -2e9, r = -2e9;
                        int derr = 0;
                        for (auto& sam : npu.waitSamples) {
                            npu.current_time = max(npu.current_time, sam.arriveTime);
                            sam.endTime = npu.current_time + sam.processTime;
                            if ( sam.userId == user_id && sam.endTime > get<1>(users[sam.userId])) {
                                    // 有超时
                                    isTimeOut = true;
                                    break;
                            }

                            if (npu.current_time > r) {
                                if (r != -2e9) {
                                    npu.schedule.emplace_back(l, r);
                                }
                                l = npu.current_time, r = sam.endTime;
                            } else {
                                r = max(r, sam.endTime);
                            }
                            
                            npu.current_time = sam.endTime;

        
                        }
                        if (r != -2e9) {
                            npu.schedule.emplace_back(l, r);
                        }      
                        // 如果有超时
                        if (isTimeOut) {
                            // 复原，并且将该用户所有段撤销并且搁置
                            npu.schedule = move(backupSchedule);
                            npu.waitSamples = move(backupWaitSample);
                            // 清空用户
                            sol.T = 0;
                            if (num >= up) {
                                timeOutCandidates.emplace_back(user_id);
                            }

                        } else {
                            // 记录用户方案
                            saveSendTimes[i] = earliestSendTime;
                            saveArriveTime[i] = arriveTime;
                            saveServers[i] = best_server + 1;
                            saveNpus[i] = best_npu + 1;
                            saveBatches[i] = batches[i];

                        }

                    }

                    if (!isTimeOut) {
                        // 接受该解
                        sol.sendTimes = move(saveSendTimes);
                        sol.arriveTimes = move(saveArriveTime);
                        sol.servers = move(saveServers);
                        sol.npus = move(saveNpus);
                        sol.batches = move(saveBatches);
                        sol.T = T;
                        solutions[user_id] = move(sol);
                        break;
                    }

                

                }
                
    
                
            }
            // 为超时的用户随机选择npu
            for (int id : timeOutCandidates) {
                int index = rng.next(0, selectOp);
                int server = allServerNpu[index].x, npu = allServerNpu[index].y;
                auto& useNpu = servers[server].npus[npu];
                // 放入到该npu等待队列的尾巴中
                int lastArriveTime  = useNpu.waitSamples.back().arriveTime;
                // 生成批次和请求
                int max_batch = min((servers[server].m - b) / a, 1000);
                int cnt_i = get<2>(users[id]);
                int s_i = get<0>(users[id]);
                int batch = min(max_batch, cnt_i);
                int T = (cnt_i + batch - 1) / batch;
                vector<int> batches(T, batch);
                batches.back() = cnt_i - (T-1)*batch;


                int startSendTime = max(s_i, lastArriveTime + 1 - latency[server][id]) ;
                int fixedLatency = latency[server][id] + 1;
                UserSolution sol;
                sol.id = id;
                sol.T = T;

                for (int i = 0;  i < T; ++ i) {
                    int earliestSendTime = startSendTime;
                    int arriveTime = earliestSendTime + latency[server][id];
                    sol.sendTimes.emplace_back(earliestSendTime);
                    sol.arriveTimes.emplace_back(arriveTime);
                    sol.servers.emplace_back(server + 1);
                    sol.npus.emplace_back(npu + 1);
                    sol.batches.emplace_back(batches[i]);
                    
                    startSendTime =  earliestSendTime + fixedLatency;
                }

                solutions[id] = sol;
            }


            function<void(const vector<UserSolution>&)> constructQueue = [&](const vector<UserSolution>& solutions)->void {
                // 根据解决方案构造队列
                retSetAllNpuState(true);
                int cnt1 = 0, cnt2 = 0;
                for (const auto& user : solutions) {
                    if (user.T != 0) ++ cnt1;
                    cnt2 ++;
                    int sz = user.batches.size();

                    assert(user.batches.size() == user.T && user.servers.size() == user.T && user.sendTimes.size() ==user.T && user.npus.size() == user.T);
                    for (int i = 0; i < sz; ++ i) {
                        int server = user.servers[i] - 1;
                        int npu = user.npus[i] - 1;
                        Request req{user.sendTimes[i] + latency[server][user.id], expense[server][user.batches[i]], 0, user.id, user.batches[i], i};
                        servers[server].npus[npu].pending.emplace(req);
  
                    }
                }
            };


            constructQueue(solutions);
            simulateReal(solutions);
            double curScore = getScore(solutions, timeOutCnt);
            if (curScore > secondBestScore) {
                secondBestScore = curScore;
                secondBestSolutions = solutions;

            }
            
            // 随机插入算子
            function<void(void)> insertNeigh = [&](void) -> void {
                // 随机选择用户
                int userId = rng.next(0, M);
                // 随机选择服务器
                int oldServer = solutions[userId].servers.size() ? solutions[userId].servers[0] - 1 : -1;
                int oldNpu = oldServer != -1 && solutions[userId].npus.size() ? solutions[userId].npus[0] - 1 : -1;
                if (oldServer == -1 || oldNpu == -1) return;

                int serverId = rng.next(0, N);
                int npu = rng.next(0, servers[serverId].g);
                while (serverId == oldServer && npu == oldNpu && runtime() < LIMIT) {
                    serverId = rng.next(0, N);
                    npu = rng.next(0, servers[serverId].g);
                }
                if (serverId == oldServer && npu == oldNpu) return;



                auto& sol = solutions[userId];

                // 保留旧数据
                int oldT = sol.T;
                vector<int> oldBatches = move(sol.batches);
                vector<int> oldServers = move(sol.servers);
                vector<int> oldNpus = move(sol.npus);
                vector<int> oldSendTimes = move(sol.sendTimes);



                // 找到新的服务器和npu，更新每个选择用户的Solution
                // 生成批次和请求
                // 优化后（高效）

                int max_batch = min((servers[serverId].m - b)/a, 1000);
                int s_i = get<0>(users[userId]);
                int cnt_i = get<2>(users[userId]);
                int batch = min(max_batch, cnt_i);

                int T = (cnt_i + batch - 1) / batch;
                sol.T = T;
                vector<int> tmp_batches(T, batch);     // 一次构造
                tmp_batches.back() = cnt_i - (T-1)*batch;       // 原地修改
                sol.batches = move(tmp_batches);  // O(1) 移动
                sol.servers = move(vector<int>(T, serverId + 1));
                sol.npus = move(vector<int>(T, npu + 1));

                // 修改发送时间
                sol.sendTimes.reserve(T);
                const int step = latency[serverId][userId] + 1;
                int current = s_i;
                int otherCost = 0;
                for (int i = 0; i < T; ++ i) {
                    sol.sendTimes.emplace_back(current);
                    current += step;
                }
                // 构造请求到队列中
                constructQueue(solutions);
                simulateReal(solutions);
                double curScore = getScore(solutions, timeOutCnt);
                // cerr << curScore << endl;
                if (curScore > secondBestScore) {
                    secondBestScore = curScore;
                    secondBestSolutions = solutions;

                } else {
                    // 撤销
                    sol.T = oldT;
                    sol.sendTimes = move(oldSendTimes);
                    sol.servers = move(oldServers);
                    sol.npus = move(oldNpus);
                    sol.batches = move(oldBatches);
                }

            };
            // 交换领域算子
            function<void(void)> swapNeigh = [&](void) -> void {
                // 随机选择两个用户
                int userId1 = rng.next(0, M);
                int userId2 = rng.next(0, M);

                if (userId1 == userId2) return;
                int server1 = solutions[userId1].servers.size() ? solutions[userId1].servers[0] - 1 : -1;
                int npu1 = solutions[userId1].npus.size() ? solutions[userId1].npus[0] - 1 : -1;
                int server2 = solutions[userId2].servers.size() ? solutions[userId2].servers[0] - 1 : -1;
                int npu2 = solutions[userId2].servers.size() ? solutions[userId2].npus[0] - 1 : -1;

                if (server1 == server2 && npu1 == npu2) return;

                vector<int> userGroups({userId1, userId2});
                vector<PII> ServerNpuGroups({{server1, npu1}, {server2, npu2}});


                // 交换数据（服务器和npu)(userId1 : (server2, npu2)); (userId2 : (server1, npu1))
                for (int i = 0; i < 2; ++ i) {
                    int another = i ^ 1;
                    auto& sol =  solutions[userGroups[i]];
                    int serverId = ServerNpuGroups[another].x, npu = ServerNpuGroups[another].y;
                    int max_batch = min((servers[serverId].m - b)/a, 1000);
                    int s_i = get<0>(users[userGroups[i]]);
                    int cnt_i = get<2>(users[userGroups[i]]);
                    int batch = min(max_batch, cnt_i);

                    int T = (cnt_i + batch - 1) / batch;
                    sol.T = T;
                    vector<int> tmp_batches(T, batch);     // 一次构造
                    tmp_batches.back() = cnt_i - (T-1)*batch;       // 原地修改
                    sol.batches = move(tmp_batches);  // O(1) 移动
                    sol.servers = move(vector<int>(T, serverId + 1));
                    sol.npus = move(vector<int>(T, npu + 1));

                    // 修改发送时间
                    sol.sendTimes.resize(T);
                    const int step = latency[serverId][userGroups[i]] + 1;
                    int current = s_i;
                    for (int i1 = 0; i1 < T; ++ i1) {
                        sol.sendTimes[i1] = (current);
                        current += step;
                    }

                }

                // 构造请求到队列中
                constructQueue(solutions);
                simulateReal(solutions);
                double curScore = getScore(solutions, timeOutCnt);
                if (curScore > secondBestScore) {
                    secondBestScore = curScore;
                    secondBestSolutions = solutions;

                } else {
                    // 撤销
                    for (int i = 0; i < 2; ++ i) {
                        auto& sol =  solutions[userGroups[i]];
                        int serverId = ServerNpuGroups[i].x, npu = ServerNpuGroups[i].y;
                        int max_batch = min((servers[serverId].m - b)/a, 1000);
                        int s_i = get<0>(users[userGroups[i]]);
                        int cnt_i = get<2>(users[userGroups[i]]);
                        int batch = min(max_batch, cnt_i);

                        int T = (cnt_i + batch - 1) / batch;
                        sol.T = T;
                        vector<int> tmp_batches(T, batch);     // 一次构造
                        tmp_batches.back() = cnt_i - (T-1)*batch;       // 原地修改
                        sol.batches = move(tmp_batches);  // O(1) 移动
                        sol.servers = move(vector<int>(T, serverId + 1));
                        sol.npus = move(vector<int>(T, npu + 1));

                        // 修改发送时间
                        sol.sendTimes.resize(T);
                        const int step = latency[serverId][userGroups[i]] + 1;
                        int current = s_i;
                        for (int i1 = 0; i1 < T; ++ i1) {
                            sol.sendTimes[i1] = (current);
                            current += step;
                        }

                    }
                }

            };

            unordered_set<int> flags;
            // 超时推理请求尝试往前移动
            function<void (vector<UserSolution> )> dealTimeOut = [&](vector<UserSolution> solutions) -> void {
                // 构造与当前解有关队列(为了将真实请求落到npu上)
                constructQueue(solutions);
                simulateReal(solutions);
                vector<UserSolution> temp = solutions;
                int processIndex = 0;
                // 遍历每个npu处理可能超时
                for (int s = 0; s < N; ++ s) {
                    for (int n = 0; n < servers[s].g; ++ n) {
                        auto &npu = servers[s].npus[n];
                        int sz = npu.schedule.size();

                        for (int j = 0; j < sz; ++ j) {
                            // 得到当前时间段的信息
                            auto [userId, batchIndex, arriveTime] = npu.infos[j];
                            auto& sol = solutions[userId];
                            int sendTime = sol.sendTimes[batchIndex];
                            int endTime = sol.endTimes[batchIndex];
                            int latestTime = get<1>(users[userId]);

                            if (endTime <= latestTime) continue;
  
                            // 当前段已经超时，计算最早到达时刻
                            int earliestSendTime = -1;
                            int preSendTime = (batchIndex > 0 ? sol.sendTimes[batchIndex - 1] : -1);
                            if (preSendTime == -1) {
                                // 当前段是第一个
                                earliestSendTime = get<0>(users[userId]);
                            } else {
                                int preServer = sol.servers[batchIndex - 1] - 1;
                                earliestSendTime = preSendTime + latency[preServer][userId] + 1;
                            }
                            // 计算最早到达时间
                            int curServer = sol.servers[batchIndex] - 1;
                            int earliestArriveTime = earliestSendTime + latency[curServer][userId];

                            if (earliestArriveTime + expense[curServer][sol.batches[batchIndex]] > latestTime ) continue;
       
                            // 往前遍历
                            int insertIndex = -1;
                            for (int j1 = j; j1 >= 0; -- j1) {
                                // 选择第一个调度区间可以放入该段，并且不超时
                                auto [l, r] = npu.schedule[j1];
                                if (l + expense[curServer][sol.batches[batchIndex]] <= latestTime) {
                                    // 找到了可以插入的区间
                                    insertIndex = j1;
                                    break;
                                }
                            }
                            if (insertIndex == -1) continue;
                            // [insertIndex, j]
                            int curIndex = j;
                            // 得到占有该区间的段信息
                            auto [insertUserId, insertBatchIndex, insertArriveTime] = npu.infos[insertIndex];
                            // 我需要调整earListArriveTime < insertArriveTime
                            // 所以依次遍历从inertIndex到curIndex中的所有用户段
                            // 调整所有到达时间大于earListArriveTime的用户段
                            // 创建副本保存原先的状态
                            unordered_map<int, vector<int>> sendTimeGroups;
                        

                            // 调整自己在最早时刻发送
                            solutions[userId].sendTimes[batchIndex] = earliestSendTime;
                            for (int j1 = insertIndex; j1 < curIndex; ++ j1) {
                                auto [changeUserId, changeBatchIndex, changeArriveTime] = npu.infos[j1];
                                sendTimeGroups[changeUserId] = (solutions[changeUserId].sendTimes);

                                if (changeArriveTime <= earliestArriveTime) {
                                    // 往后移动
                                    // 将用户的大于等于changeBatchIndex的batch全部往后移动
                                    int needArriveTime = earliestArriveTime + 1;
                                    int needSendTime = needArriveTime - latency[s][changeUserId];
                                    needSendTime = max(needSendTime, solutions[changeUserId].sendTimes[changeBatchIndex]);

                                    // 对该用户的后续段进行挑选
                                    int move = needSendTime - solutions[changeUserId].sendTimes[changeBatchIndex];
                                    solutions[changeUserId].sendTimes[changeBatchIndex] = needSendTime;

                                    for (int j2 = changeBatchIndex + 1; j2 < solutions[changeUserId].T; ++ j2) {
                                        solutions[changeUserId].sendTimes[j2] = max(solutions[changeUserId].sendTimes[j2], solutions[changeUserId].sendTimes[j2 - 1] + latency[s][changeUserId] + 1);

                                    }


                                }
                            }

                            // 处理完一个超时计算一下请求
                            constructQueue(solutions);
                            simulateReal(solutions);
                            double score = getScore(solutions, timeOutCnt);
                            if (score > secondBestScore) {
                                secondBestScore = score;
                                secondBestSolutions = solutions;

                                return;

                            } else {
                                
                                solutions = temp;
                                constructQueue(solutions);
                                simulateReal(solutions);
                            }


                        }

                    }
                }
                
            } ;
            // dealTimeOut(solutions);
            // 随机交换服务器和npu
            // startTime - LIMIT : startTime + (LIMIT - startTime) / 2
            int useLimit = LIMIT;
            if (is == 0) useLimit = startTime + (LIMIT - startTime) / 2;
            function<void(int)> climb = [&](int iterTime) -> void{
                while (runtime() < useLimit) {
                    swapNeigh();
                    insertNeigh();
                    // solutions = secondBestSolutions;
                    // dealTimeOut(solutions);


                }
            };

            climb(LIMIT);


        }
        int timeOutCnt = 0;
        double score = getScore(secondBestSolutions, timeOutCnt);
        if (score > bestScore) {
            bestScore = score;
            bestSolutions = secondBestSolutions;
        }

        // cerr << "最好的分数: " << bestScore << endl;


    }

    
    void mainLoop() {
        int maxSegment = 300;
        // 计算每个用户对于所有服务器的最小分段数
        vector<vector<int>> minSegment(N, vector<int>(M, 1));
        for (int s = 0; s < N; ++ s) {
            for (int u = 0; u < M; ++ u) {
                int maxBatch = min((servers[s].m - b) / a, 1000);
                int cnt = get<2>(users[u]);
                // 遍历分段
                int segment = 1;
                while (true) {
                    if (maxSegment * (maxBatch / segment) > cnt) {
                        ++ segment;
                    } else {
                        -- segment;
                        break;
                    }
                }
                minSegment[s][u] = max(segment, minSegment[s][u]);
            }
        }
        // 方案一(枚举分段)
        for (int segment = 2; segment <= 9; ++ segment) schedulingOne(minSegment, segment);

        int secondStartTime = runtime();
        // 方案二
        schedulingTwo(secondStartTime);


    }
    // 计算处理x个任务所需的时间
    int computeTime(int x, int k) {
        if (x <= 0) return 0;
        double time = sqrt((double)x / k);
        // 向上取整
        return (int)ceil(time);
    }

    // 找到最优的任务分配序列并返回最小总时间
    int findOptimalSchedule(int cnt, int k, int maxCnt, vector<int>& schedule) {
        // dp[i]表示处理前i个任务的最小时间
        vector<int> dp(cnt + 1, INT_MAX);
        // prev[i]记录处理前i个任务的最后一步处理了多少个任务
        vector<int> prev(cnt + 1, 0);
        
        dp[0] = 0; // 处理0个任务时间为0
        
        for (int i = 1; i <= cnt; ++i) {
            // 尝试最后一步处理j个任务，j的范围是1到min(maxCnt, i)
            for (int j = 1; j <= min(maxCnt, i); ++j) {
                int prevCount = i - j;
                if (dp[prevCount] != INT_MAX) {
                    int currentTime = dp[prevCount] + computeTime(j, k);
                    if (currentTime < dp[i]) {
                        dp[i] = currentTime;
                        prev[i] = j;
                    }
                }
            }
        }
        
        // 回溯构建任务分配序列
        int remaining = cnt;
        while (remaining > 0) {
            schedule.push_back(prev[remaining]);
            remaining -= prev[remaining];
        }
        reverse(schedule.begin(), schedule.end());
        
        return dp[cnt];
    }
    INLINE double getScore(const vector<UserSolution> &userSolution, int& timeOutCnt){
        timeOutCnt = 0;
        double score = 10000.0;
        int sz = userSolution.size();
        int cnt = 0;
        for (const auto& user : userSolution) {
            cnt += (user.endTime > std::get<1>(users[user.id]));
            #ifdef TIMEOUTINFO
                if ((user.endTime > std::get<1>(users[user.id]))) {
                    cerr << "********************" << endl;
                    cerr << "超时用户: " << user.id << " endTime: " << user.endTime <<  endl;
                    cerr << get<0>(users[user.id]) << " " << get<1>(users[user.id]) << endl;
                    for (int i = 0; i < user.batches.size(); ++ i) {
                        cerr << "sendTime: " << user.sendTimes[i] << " endTime: " << user.endTimes[i] << " server: " << user.servers[i] - 1<< " npu: " << user.npus[i] - 1<< 
                        " runTime: " << expense[user.servers[i] - 1][user.batches[i]] << endl;
                    }
                }
            #endif

        }
        timeOutCnt = cnt;
        score *= hxInt[cnt];
        double sum = 0.0;
        for(const auto& user : userSolution) {
            int s = std::get<0>(users[user.id]), e = std::get<1>(users[user.id]);
            double item1 = h(((user.endTime - e) / (1.0 * e - s)));
         
            double item2 = pxInt[user.moveCount];
            // 遍历
            sum += item1 * item2;

        }
        
        score *= sum;


        return score;


    }
    // 真实调度模拟
    void simulateReal(vector<UserSolution> &userSolution) {
        // 状态
        retSetAllNpuState(false);
        // 重置用户的结束时间
        for (auto& user : userSolution) user.endTime = 0, user.valid = false, user.endTimes.assign(user.batches.size(), 0);
        
        int totalBatch = 0;
        // 第二阶段：处理所有NPU的请求队列
        for (int s = 0; s < N; ++ s) {
            for (int n = 0; n < servers[s].g; ++ n) {
                auto& npu = servers[s].npus[n];
                while (!npu.pending.empty()) {
                    Request req = npu.pending.top();
                    npu.pending.pop();
                    
                    // 时间推进到请求到达时间
                    if (req.arrive_time > npu.current_time) {
                        npu.current_time = req.arrive_time;
                    }          
                    // 处理请求
                    int end_time = npu.current_time + req.process_time;
                    
                    userSolution[req.user_id].endTimes[req.index] = end_time;
                    
                    // 记录时间窗口
                    npu.schedule.emplace_back(npu.current_time, end_time);
                    npu.infos.emplace_back(req.user_id, req.index, req.arrive_time);
                    // 更新用户结束时间
                    userSolution[req.user_id].endTime = max(
                        userSolution[req.user_id].endTime, end_time);
                    
                    // 释放显存（实际应记录释放时间）
                    npu.current_time = end_time;

                }
            }
        }

        for (auto& sol : userSolution) {
            int e_i = get<1>(users[sol.id]), cnt_i = get<2>(users[sol.id]);
            sol.valid = (accumulate(sol.batches.begin(), sol.batches.end(), 0) == cnt_i);

            sol.valid &= all_of(sol.endTimes.begin(), sol.endTimes.end(), [&](int t){
                return t != 0;
            });

        }

        #ifdef NPUINFO
            // 输出每个npu的调取区间
            for (int s = 0; s < N; ++ s) {
                cerr << "server: " << s << endl;
                for (int n = 0; n < servers[s].g; ++ n) {
                    cerr << "npu: " << n << endl;
                    auto& npu = servers[s].npus[n];
                    for (const auto& [l, r] : npu.schedule) {
                        cerr << "[" << l << " " << r << "]";
                    }
                    cerr << endl;
                }
            }
        #endif
    }
    INLINE double calculateUpperScore() {
        double score = 10000;
        return 1e4 * (M * (1.0 * pow(2, 1.0 / 100)));
    }

    INLINE void retSetAllNpuState(bool flag) {
        for (int s = 0; s < N; ++ s) {
            for (int n = 0; n < servers[s].g; ++ n) {
                auto& npu = servers[s].npus[n];
                npu.reset(flag);
            }
        }
    }

    void output() {
        for (int i = 0; i < M; ++i) {
            auto& sol = bestSolutions[i];
            printf("%d\n", sol.T);
            for (size_t j = 0; j < sol.T; ++j) {
                printf("%d %d %d %d", 
                    sol.sendTimes[j], 
                    sol.servers[j], 
                    sol.npus[j], 
                    sol.batches[j]);
                if (j != sol.T-1) putchar(' ');
            }
            putchar('\n');
        }
    }
};

int main() {


    #ifdef LOCAL
        if (fopen("../data.in", "r") != nullptr) {

            freopen(("../input.txt"), "r", stdin);
            freopen(("./output.txt"), "w", stdout);

            
        } else {
            cerr << "open file uncorrectly" << endl;
        }
    #endif 
    preProcess();

    static Scheduler scheduler;
    scheduler.init();
    scheduler.mainLoop();
    scheduler.output();
    return 0;
}