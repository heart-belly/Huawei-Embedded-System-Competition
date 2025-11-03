
#include <iostream>
#include <cstring>
#include <algorithm>
#include <vector>
#include <queue>
#include <unordered_set>
#include <set>
#include <random>
using namespace std;



// 用于改变通道寻路
struct EdgeInfo {
    double distance;
    int preIndex, ver, edgeId, channelL, time;
    EdgeInfo(double _distance, int _preIndex, int _ver ,int _edgeId, int _channelL, int _time) {
        distance = _distance;
        preIndex = _preIndex;
        ver = _ver;
        edgeId = _edgeId;
        channelL = _channelL;
        time = _time;
    }
};
typedef pair<double, EdgeInfo> PDE;
typedef pair<int, int> PII;
typedef pair<int, pair<pair<int, int>, pair<int, int>>> PIPII;
typedef pair<pair<int, int>, pair<int, int>> PIIPIIP;
typedef pair<double, int> PDI;

#define x first
#define y second

const int N = 210, M = 1010, K = 40; // K : 边上的通道数量
const int JJ = 5010;


struct cmp
{
	bool operator()(const pair<double, EdgeInfo>& node1, const pair<double, EdgeInfo>& node2)
	{
        if (node1.x != node2.x ) {
            return node1.x > node2.x;
        }

        return node1.y.time > node2.y.time;

	}
};
/*
* 场景的初始状态
*/

unordered_set<int> _edgeHasBusiness[M]; // 当前边已有的业务
int _edgeChannel[M][K + 1]; // 边上的通道的使用情况

// 节点变通数相关存储
struct P {
    int curUseP, remainUseP; // curUseP : 当前已用变通数; remainUseP : 剩余可用变通数; sumP : curUseP + remainUseP;
} p[N], _p[N];
// 存储节点对应的包含变通道的业务
vector<unordered_set<int>> nodeChangeHasWork;
// 业务相关存储
struct BUSINESS {
    int src; // 起点
    int snk; // 终点
    int s; // 经过的边数
    int width; // 通道的宽度
    int v; // 该业务的价值
    int l, r; // 该业务起点出发占有的通道范围
    int useChange; // 使用的变通道数量

    vector<vector<int>> ownEdges; // 0 ： 边的编号; 1 : 通道的开始编号; 2 : 通道的结束编号
} business[JJ], _business[JJ];
/*
边的通道数的实时变化
边的顺序问题
*/
vector<int> g[N][N]; // 构造邻接矩�?(以点为下标作为索引边的集�? -- 根据输入顺序构造边的存入顺�?) *** 由于是无向边，需要保证数据的一致�?
PII edgeToNode[M]; // 根据边的编号找到是在哪两个结点之�? (注意的出来的结点不分先后，因为是无向图)
unordered_set<int> edgeHasBusiness[M]; // 当前边已有的业务
bool edgeHasBusinessReal[M][JJ];
bool nodeChangeHasWorkReal[N][JJ];
int edgeChannel[M][K + 1]; // 边上的通道的使用情况

/*
* 寻求点与点之间的路径
*/
vector<vector<int>> nodeToNode;
vector<vector<int>> nodeToEdge;
int n, m; // n : 节点; m : 边数;
int J; // 表示业务数量 (1 <= J <= 5000)
int T; // 表示测试场景

/*
* 寻找最短路径
*/
bool st[M]; // 标记边有没有被遍历过
bool stt[M][K + 1]; // 边 + 通道起点作为唯一标识
int dist[N], cnt[N];
PIIPIIP pre3[M]; // 1 : preVer ; 2 : preEdge ; 3 : l ; 4 : r (用来构造路径)
PIIPIIP pre2[M * 40]; // 保存（边 + 通道）的信息 {1 ：起点; 2 : 终点; 3 : 边的编号; 4 : 通道l}
int pre[M * 40]; // 保存上一条（边 + 通道）的唯一标识
int idx, ttime;
PII canUseChannel[100]; // x : 通道起点; y : 赋予该边和通道唯一的权值
int idxCanUse;
unordered_set<int> memo[M];
bool memoSt[M][K + 1];


// 保存已经坏边以及死亡的业务
unordered_set<int> badEages, badBusiness;
// 回答数组
vector<PII> answer;
unordered_set<int> succeedWork;

/*
* 资源分配数据结构
*/

vector<vector<int>> clearEdgeInfoVec; // workId, e, l, r
vector<vector<int>> changeEdgeInfoVec; // workId, e, l, r
vector<PII> subNodeAlternate; // 减少节点的变通道次数
vector<PII> addNodeAlternate; // 增加节点的变通道次数

vector<set<int>> saveBadEdge;

int ptotal = 0;
int total = 0, totalValue = 0;
int curDamage = 0;

int T1;
/*
*/
/*
* 初始化
*/
void init(); // 程序初始化
void input(); // 输入有关的全局信息
void initData(); // 初始化每次场景的数据
void interation(); // 进行交互的操作
void initFindPathDaTa(int workId); // 初始化寻路时的相关数据
/*
* 功能性函数
*/
void dijkstra(int end); // 求单源最短路径
int getCanUseCannel(int edgeId);
PII selectChannel(bool isChange, int width, int prel, int prer, int edgeId); // 选择边上的通道(1:该边是否可以使用变通道; 2:通道的宽度; (3,4) : 该边前一条边的通道使用情况)
void changeReInfo(int workId, int edgeId, int l, int r); // 保持数据的一致性
void clearEdgeInfo(int workId, int edgeId, int l, int r); // 保持数据一致性
double jaccardSimilarity(const set<int>& set1, const set<int>& set2); // 计算两个集合的Jaccard相似度系数
set<int> generateSubset(const vector<int>& originalSet, int subsetSize); // 生成一个随机的子集
void dataProcessing(); // 处理边的划分集合
void constructTest1(); // 构造测试集1
void constructionPre3(int index, int workId); // 构造路径数组
std::vector<std::set<int>> heuristicSearch(const std::vector<std::set<int>>& sets, const std::set<int>& originalSet, int subsetSize, int numSubsets);

/*
* 搜索路径时，针对不同需求设计不同的策略
*/
bool satisfyFix(int edgeId, int workId, int preNode, int nextNode, int preL, int preR); // 判断当前边能不能满足业务需求
bool satisfyChange(int &index, int &workId, int &edgeId); // 判断当前边能不能满足业务需求

/*
* 处理策略
*/
void deal1(int start, int end, int workId, int edgeId); // 寻找受损边两端不变通道的路径
int aStar1(int wordId, int start, int end, int preL, int preR);
void deal2(int start, int end, int workId, int op); // 寻找起点到终点不变通道的路径
int aStar2(int wordId, int start, int end, int preL, int preR);
void deal3(int start, int end, int workId, int op); // 寻找起点到终点可变通道的路径(k : 表示传入变通道的权重)
int aStar3(int wordId, int start, int end);

/*
* 资源分配相关函数
*/
void changeReInfo(int workId, int edgeId, int l, int r); // 保持数据的一致性
void clearEdgeInfo(int workId, int edgeId, int l, int r); // 保持数据一致性
void dealEdgeAndNode();

int main() {
    init();
    input();
    initData();
    dataProcessing();
    constructTest1();
    // cout << 0 << endl;
    // fflush(stdout);
    cin >> T;
    for (int t = 0; t < T; t ++ ) {
        // 初始化数据
        initData();
        if (t < T1) {
            for (int mayEdgeId : saveBadEdge[t]) {
                badEages.insert(mayEdgeId);
            }
        }
        // 进入交互部分
        interation();

        for (int i = 1; i <= J; i ++ ) {
            if (!badBusiness.count(i)) {
                totalValue += business[i].v;
                total += 1;
            }
        }
    }


    cout << "total value :" << totalValue << endl;
    cout << "total :" << total << endl;
    cout << "score : " << (int)(10000.0 * totalValue / ptotal) << endl;


}

void init() {
    // 关闭同步输入流
    std::ios::sync_with_stdio(false);
    cin.tie(0);
    cout.tie(0);

}
void initData() {
    // 将初始化场景赋值给使用场景
    for (int i = 1; i <= n; i ++ ) {
        p[i].curUseP = _p[i].curUseP;
        p[i].remainUseP = _p[i].remainUseP;
    }
    for (int i = 1; i <= m; i ++ ) {
        for (int j = 1; j <= K; j ++ ) {
            edgeChannel[i][j] = _edgeChannel[i][j];
        }
    }
    // 业务赋值
    for (int i = 1; i <= J; i ++ ) {
        business[i] = _business[i];
    }
    // 边和业务的赋值
    memset(edgeHasBusinessReal, 0, sizeof edgeHasBusinessReal);
    memset(nodeChangeHasWorkReal, 0, sizeof nodeChangeHasWorkReal);
    for (int i = 1; i <= m; i ++ ) {
        edgeHasBusiness[i] = _edgeHasBusiness[i];
        for (auto work : edgeHasBusiness[i]) {
            edgeHasBusinessReal[i][work] = true;
        }
    }

    // 坏边的集合刷新
    badEages = {};
    // 坏的业务刷新
    badBusiness = {};

    // 刷新通过节点改变通道的业务
    nodeChangeHasWork = vector<unordered_set<int>> (n + 1);
}
void input() {
    cin >> n >> m; // 输入结点和边
    // 输入结点的最大变通数
    for (int i = 1; i <= n; i ++ ) {
        cin >> _p[i].remainUseP;
        // 完成初始化
        _p[i].curUseP = 0;


    }
    nodeToNode = vector<vector<int>> (n + 1);
    nodeToEdge = vector<vector<int>> (n + 1);
    // 考虑重边情况（节点之间的重边是可以使用的
    for (int i = 1; i <= m; i ++ ) {
        int a, b;
        cin >> a >> b;
        // 以在数组中的下标作为边的唯一索引
        g[a][b].push_back(i);
        g[b][a].push_back(i);

        // 保存点可以到达的点
        nodeToNode[a].push_back(b);
        nodeToNode[b].push_back(a);

        nodeToEdge[a].push_back(i);
        nodeToEdge[b].push_back(i);

        edgeToNode[i] = {a, b};
        for (int j = 1; j <= K; j ++ ) {
            // -1 : 表示可用通道
            _edgeChannel[i][j] = -1;
        }
    }
    // 输入图中的业务(对边的处理以及业务的处理)
    cin >> J;
    for (int i = 1; i <= J; i ++ ) {
        int l, r;
        cin >> _business[i].src >> _business[i].snk >> _business[i].s >> l >> r >> _business[i].v;
        _business[i].width = (r - l + 1);
        _business[i].l = l, _business[i].r = r;
        _business[i].useChange = 0;
        ptotal += _business[i].v;
        vector<vector<int>> temp;
        for (int j = 0; j < _business[i].s; j ++ ) {
            int e;
            cin >> e;
            _edgeHasBusiness[e].insert(i); // e边当前有业务i在上面运行
            auto [a, b] = edgeToNode[e];

            temp.push_back({e, l, r});

            // 遍历边的通道
            for (int j1 = l; j1 <= r; j1 ++ ) {
                // 通道被业务i占领
                _edgeChannel[e][j1] = i;
            }

        }

        _business[i].ownEdges = temp;
    }
    for (int i = 1; i <= n; i ++ ) {
        sort(nodeToEdge[i].begin(), nodeToEdge[i].end(), [&](const int &a, const int &b){
            double rate1 = 0.0, rate2 = 0.0;
            for (int k1 = 1; k1 <= K; k1 ++ ) {
                if (edgeChannel[a][k1] != -1 ) {
                    rate1 += 1.0 / 40;
                }
                if (edgeChannel[b][k1] != -1) {
                    rate2 += 1.0 / 40;
                }
            }

  
            return rate1 < rate2;
        });
    }

}

void interation() {
    int op; // 坏的边
    while (cin >> op, op != -1) {
        // 找到在该边上运行的业务数
        auto &workVec = edgeHasBusiness[op];
        curDamage = 0;
        for (int id : workVec) {
            if (badBusiness.count(id)) curDamage ++;
        }
        //将坏的边加入到集合中
        badEages.insert(op);
        // 刷新一下回答数组
        answer.clear();
        succeedWork = {};
        clearEdgeInfoVec = {}, changeEdgeInfoVec = {};
        subNodeAlternate = {}, addNodeAlternate = {};
        // 寻找下在改边上运行的业务，保证业务必须活着
        // 该边上没有业务在运行
        if (workVec.empty()) {
            cout << 0 << endl;
            fflush(stdout);
        } else {
            // 找到在该边上运行的业务根据价值排序
            vector<pair<double, int>> sortG;
            for (auto id : workVec) {
                int start = business[id].src;
                int end = business[id].snk;
                sortG.push_back({pow(1.0 * business[id].v, 1) * pow(1.0 * business[id].ownEdges.size(), 0.4) * pow(1.0 * business[id].width, 0.6) * pow(1.0 * max(1, business[id].useChange), 0.1), id});

            }

            sort(sortG.begin(), sortG.end(), [&] (const PDI& a, const PDI& b) {
                return a.x > b.x;
            });
            // 处理不改变重边所对应的通道
            for (int i = 1; i <= 2; i ++ ) {
                for (auto [_, workId] :  sortG) {
                    if (!edgeHasBusinessReal[op][workId]) continue;
                    // 当前业务已死不用规划或者前一次已成功规划
                    if (badBusiness.count(workId) || succeedWork.count(workId)) continue;

                    // 找到业务的起点和终点
                    int start = business[workId].src, end = business[workId].snk;

                    if (i == 1) {
                        deal3(start, end, workId, op);
                    }  else if (i == 2) {
                        badBusiness.insert(workId);  
                    } 

                }

            }
            // 统一处理边与节点的资源
            // cout << "have fluenced :" << workVec.size() - curDamage << endl;
            // 统一处理边与节点的资源
            dealEdgeAndNode();
            // cout << "survive :";
            cout << answer.size() << endl;
            fflush(stdout);
            for (int i = 0; i < answer.size(); i ++ ) {
                // first : 业务编号; second : 经过的边的数量
                // cout << "my answer : ";
                // cout << "sdsad" << business[answer[i].x].width << endl;
                cout << answer[i].x << " " << answer[i].y << endl;
                fflush(stdout);
                auto backup = business[answer[i].x].ownEdges;
                for (int j = 0; j < backup.size(); j ++ ) {
                    cout << backup[j][0] << " " << backup[j][1] << " " << backup[j][2] << " ";
                }
                cout << endl;
                fflush(stdout);
            }

            // 该边已经死了，清空在该边上运行的业务
            workVec.clear();
        }
    }
}



void dijkstra(int end) {
    priority_queue<PII, vector<PII>, greater<PII>> heap;
    memset(st, 0, sizeof st);
    memset(dist, 0x3f, sizeof dist);
    dist[end] = 0;
    heap.push({0, end});
    
    while (heap.size()) {
        auto t = heap.top();
        heap.pop();
        
        int ver = t.y, dis = t.x;
        
        if (st[ver]) continue;
        
        st[ver] = true;

        for (int i = 1; i <= n; i ++ ) {
            if (!st[i]) {
                // 判断有没有变可以走
                bool flag = false;
                for (int edgeId : g[ver][i]) {
                    if (!badEages.count(edgeId)) {
                        flag = true;
                        break;
                    }
                }

                if (flag && dist[i] > dist[ver] + 1) {
                    dist[i] = dist[ver] + 1;
                    heap.push({dist[i], i});
                }
            }
        }
    }
}




bool satisfyChange(int &index, int &workId, int &edgeId) {
    // 边如果之前以及坏了就不能走
    if (badEages.count(edgeId)) return false;
    // 得到业务相关信息
    int width = business[workId].width;

    if (memo[edgeId].size()) {
        int preNode = pre2[index - 1].x.x, nextNode = pre2[index - 1].x.y;
        int  channel = pre2[index - 1].y.y, preEdgeId = pre2[index - 1].y.x;
        
        if (memo[edgeId].count(channel)) {
            canUseChannel[idxCanUse ++ ] = {channel, 1};
            memoSt[edgeId][channel] = true;
        }

        for (int ll : memo[edgeId]) {
            if (ll != channel && !memoSt[edgeId][ll] && (nodeChangeHasWorkReal[nextNode][workId] || p[nextNode].remainUseP > 0)) {
                    if (nodeChangeHasWorkReal[nextNode][workId]) {
                        canUseChannel[idxCanUse ++] = {ll, 3};
                    } else {
                        canUseChannel[idxCanUse ++] = {ll, 3 + 5 / p[nextNode].remainUseP};
                    }
                memoSt[edgeId][ll] = true;
            }
        }

    } else {
        PII temp[50];
        int idxTemp = 0;
        int acc = 0;
        for (int r = 1; r <= 40; r ++ ) {

            if (edgeChannel[edgeId][r] == -1 || edgeChannel[edgeId][r] == workId) {
                acc ++;
            } else if (acc > 0){
                temp[idxTemp ++] = {r - acc, r - 1};
                acc = 0;
            }

        }
        // 处理边界问题
        if (edgeChannel[edgeId][40] == -1 || edgeChannel[edgeId][40] == workId) {
            temp[idxTemp ++] = {40 - acc + 1, 40};
        }

        sort(temp + 0, temp + idxTemp, [&](const PII &p1, const PII &p2){
            return p1.y - p1.x < p2.y - p2.x;
        });

        // index 表示当前的点前一条（边 + 通道）
        if (index == 1) {
            // 表示从起点出发寻找的边
            for (int i = 0; i < idxTemp; i ++ ) {
                int l = temp[i].x, r = temp[i].y;
                if (r - l + 1 < width) continue;
                int ll = l;
                while (ll + width - 1 <= r) {
                    if (ll == business[workId].l) {
                        canUseChannel[idxCanUse ++] = {ll, 0};  
                    } else {
                        canUseChannel[idxCanUse ++] = {ll, 1};   
                    }
                    memo[edgeId].insert(ll);
                    memoSt[edgeId][ll] = true;
                    ll ++;
                }
            }

        } else {
            int preNode = pre2[index - 1].x.x, nextNode = pre2[index - 1].x.y;
            int  channel = pre2[index - 1].y.y, preEdgeId = pre2[index - 1].y.x;
            // 得到上一条边和通道的相关数据
            for (int i = 0; i < idxTemp; i ++ ) {
                int l = temp[i].x, r = temp[i].y;
                if (r - l + 1 < width) continue;
                int ll = l;
                while (ll + width - 1 <= r) {
                    if (ll == channel) {
                        canUseChannel[idxCanUse ++] = {ll, 1};
                        memoSt[edgeId][ll] = true;
                    } else if ((nodeChangeHasWorkReal[nextNode][workId] || p[nextNode].remainUseP > 0) ){
                        if (nodeChangeHasWorkReal[nextNode][workId]) {
                            canUseChannel[idxCanUse ++] = {ll, 3};
                        } else {
                            canUseChannel[idxCanUse ++] = {ll, 3 + 5 / p[nextNode].remainUseP};
                        }
                        memoSt[edgeId][ll] = true;
                    }
                    memo[edgeId].insert(ll);
                    ll ++;
                }
            }

        }
    }


    return idxCanUse > 0;
}




void deal3(int start, int end, int workId, int edgeId) {

    // 寻找从业务头到尾的一条简单路径（可变通道）
    auto backup = business[workId].ownEdges;
    // 从起点出发的一个通道范围
    int l = business[workId].l, r = business[workId].r;

    // 🐕
    dijkstra(end);
    int index = aStar3(workId, start, end);
    if (index != -1 ) {
        // 构造pre3
        constructionPre3(index, workId);
        // 找到路径(遍历边的编号进行处理)
        int endNode = end;
        string pathNode = "";
        vector<int> nodeTemp; // 保存路径中的节点
        vector<vector<int>> pathEdge;
        while ( endNode != -1) {
            pathNode += to_string(endNode) + " ";
            nodeTemp.push_back(endNode);
            if (pre3[endNode].x.y != -1) {
                pathEdge.push_back({pre3[endNode].x.y, pre3[endNode].y.x, pre3[endNode].y.y});
            }
            endNode = pre3[endNode].x.x;
        }

        // reverse(pathNode.begin(), pathNode.end());
        reverse(nodeTemp.begin(), nodeTemp.end());
        reverse(pathEdge.begin(), pathEdge.end());

        // 输出边和节点的路径
        // cout << "node :";
        // for (int i = 0; i < nodeTemp.size(); i ++ ) {
        //     cout << nodeTemp[i] << " ";
        // }
        // cout << endl;
        // cout << "edge :";
        // for (int i = 0; i < pathEdge.size(); i ++ ) {
        //     cout << pathEdge[i][0] << " ";
        // }
        // cout << endl;
        // cout << "node express path : " << pathNode << endl;


        // 处理起点到终点的一整条路径
        vector<vector<int>> temp;
        auto backup = business[workId].ownEdges;
        // 释放业务所用边和结点的资源（边的通道数和结点的变通道能力）

        int curL = -1, curR = -1;
        int curNode = start;
        for (int i = 0; i < backup.size(); i ++ ) {
            int e = backup[i][0], l = backup[i][1], r = backup[i][2];

            if (curL == -1 || curR == -1) {
                curL = l, curR = r;
            } else {
                if (curL != l || curR != r) {
                    addNodeAlternate.push_back({curNode, workId});
                }

                curL = l, curR = r;
            }
            // 用异或变节点
            curNode = curNode ^ edgeToNode[e].x ^ edgeToNode[e].y;
            // 修改边上运行业务的相关数据
            // clearEdgeInfo(workId, e, l, r);
            clearEdgeInfoVec.push_back({workId, e, l, r});

        }

        curL = -1, curR = -1;
        business[workId].useChange = 0;
        // cout << "edge express path :";
        for (int i = 0; i < pathEdge.size(); i ++ ) {
            
            if (curL == -1 || curR == -1) {
                curL = pathEdge[i][1], curR = pathEdge[i][2];
            } else {
                if (curL != pathEdge[i][1] || curR != pathEdge[i][2]) {
                    business[workId].useChange ++;
                    if (nodeChangeHasWorkReal[nodeTemp[i]][workId]) {
                        subNodeAlternate.push_back({nodeTemp[i], workId});
                    } else {
                        nodeChangeHasWorkReal[nodeTemp[i]][workId] = true;
                        p[nodeTemp[i]].curUseP += 1, p[nodeTemp[i]].remainUseP -= 1;
                    }
                }
                curL = pathEdge[i][1], curR = pathEdge[i][2];
            }
            // cout << pathEdge[i] << " ";
            temp.push_back({pathEdge[i][0], pathEdge[i][1], pathEdge[i][2]});
            // 与变pathEdge[i] 对应的节点为 pathNode[i]
            // 修改边上运行业务的相关数据
            changeReInfo(workId, pathEdge[i][0], pathEdge[i][1], pathEdge[i][2]);
            changeEdgeInfoVec.push_back({workId, pathEdge[i][0], pathEdge[i][1], pathEdge[i][2]});
        }



        business[workId].ownEdges = temp;

        answer.push_back({workId, temp.size()});
        // for (int i = 0; i < business[workId].ownEdges.size(); i ++ ) {
        //     cout << business[workId].ownEdges[i][0] << " ";
        // }
        // cout << endl;


        succeedWork.insert(workId);


    } 


 }


 int aStar3(int workId, int start, int end) {

    priority_queue<PDE, vector<PDE>, cmp> heap;

    initFindPathDaTa(workId);


    auto backup = business[workId].ownEdges;


    heap.push({dist[start], {0, 0, start, -1, -1, ttime}});


    while (heap.size()) {
        auto t = heap.top();
        heap.pop();
        int distance = t.y.distance, preIndex = t.y.preIndex, ver = t.y.ver, edgeId = t.y.edgeId, channelL = t.y.channelL;

        if (ver != start) {
            if (stt[edgeId][channelL]) continue;
            stt[edgeId][channelL] = true;
            st[ver] = true;
            // 求最短路径
            pre[idx] = preIndex;
            int preNode = start;
            if (preIndex != 0) {
                preNode = pre2[preIndex].x.y;
            }
            pre2[idx] = {{preNode, ver}, {edgeId, channelL}};

            idx ++;
        }


        st[ver] = true;
        if (st[end]) return idx;

        unordered_set<int> pathHasNode;
        // 寻找路径上已经经过的节点
        int temp = idx - 1;
        while (temp != 0) {
            pathHasNode.insert(pre2[temp].x.x);
            pathHasNode.insert(pre2[temp].x.y);
            temp = pre[temp];
        }
        pathHasNode.insert(start);
        if (pathHasNode.size() >= business[workId].ownEdges.size() + 5) continue;

        // for (int i : nodeToNode[ver]) {
        //     if (pathHasNode.count(i) == 0) {
        //         for (int edgeId : g[ver][i]) {
        //             idxCanUse = 0;
        //             if (satisfyChange(idx, workId, edgeId)) {
        //                 for (int j = 0; j <  idxCanUse; j ++) {
        //                     int l = canUseChannel[j].x, k = canUseChannel[j].y; 
        //                     heap.push({distance + k + dist[i], {distance + k, idx - 1, i, edgeId, l, ++ ttime}});
        //                 }
        //             }
        //         }
        //     }            
        // }
        for (int edgeId : nodeToEdge[ver]) {
            auto [a, b] = edgeToNode[edgeId];
            int node = a ^ b ^ ver;
            if (pathHasNode.count(node) == 0) {
                idxCanUse = 0;
                if (satisfyChange(idx, workId, edgeId)) {
                    for (int j = 0; j <  idxCanUse; j ++) {
                        int l = canUseChannel[j].x, k = canUseChannel[j].y; 
                        heap.push({distance + k + dist[node], {distance + k, idx - 1, node, edgeId, l, ++ ttime}});
                    }
                }
            }
        }
        // for (int i = 1; i <= n; i ++ ) {
        //     if (pathHasNode.count(i) == 0) {
        //         for (int edgeId : g[ver][i]) {
        //             idxCanUse = 0;
        //             if (satisfyChange(idx, workId, edgeId)) {
        //                 for (int j = 0; j <  idxCanUse; j ++) {
        //                     int l = canUseChannel[j].x, k = canUseChannel[j].y; 
        //                     heap.push({distance + k + dist[i], {distance + k, idx - 1, i, edgeId, l, ++ ttime}});
        //                 }
        //             }
        //         }
        //     }
        // }
    }
    
    return -1;
}



PII selectChannel(bool isChange, int width, int prel, int prer, int edgeId) {

    if (!isChange) {

        for (int i = prel; i <= prer; i ++ ) {
            if (edgeChannel[edgeId][i] != -1) {
                return {-1, -1};
            }
        }

        return {prel, prer};
    } else {
        // 先判断能不能不选择边通道
        bool flag = true;
        for (int i = prel; i <= prer; i ++ ) {
            if (edgeChannel[edgeId][i] != -1) {
                flag = false;
                break;
            }
        }
        if (flag) {
            return {prel, prer};
        }

        // 选择改变通道
        int l = -1, r = -1;
        int cur = 0;
        for (int i = 1; i <= K; i ++ ) {
            if (edgeChannel[edgeId][i] == -1) {
                cur ++;
                if (l == -1) {
                    l = i;
                }
                r = max(r, i);
            } else {
                if (cur >= width) {
                    return {l, l + width - 1};
                }

                cur = 0;
                l = -1;

            }
        }

        if (l != -1) {
            if (cur >= width) return {l, l + width - 1};
            else {
                return {-1, -1};
            }
        }
    }

    return {-1, -1};
}


void changeReInfo(int workId, int edgeId, int l, int r) {

    // 改变边的情况
    // 1 : 边上运行业务的情况
    edgeHasBusiness[edgeId].insert(workId);
    edgeHasBusinessReal[edgeId][workId] = true;
    // 2 : 边上的通道数
    for (int i = l; i <= r; i ++ ) {
        edgeChannel[edgeId][i] = workId;
    }
    // cout << "workId :" << workId << endl;
    // cout << edgeId << " " << edgeHasBusiness[edgeId].size() << endl;
}

void clearEdgeInfo(int workId, int edgeId, int l, int r) {
    // 1 : 边上运行业务的情况

    edgeHasBusinessReal[edgeId][workId] = false;
    edgeHasBusiness[edgeId].erase(workId);
     // 2 : 边上的通道数
    for (int i = l; i <= r; i ++ ) {
        edgeChannel[edgeId][i] = -1;
    }
    
}

void dealEdgeAndNode() {
    // 释放资源
    for (int i = 0; i < clearEdgeInfoVec.size(); i ++ ) {
        int workId = clearEdgeInfoVec[i][0], edgeId = clearEdgeInfoVec[i][1], l = clearEdgeInfoVec[i][2], r = clearEdgeInfoVec[i][3];
        clearEdgeInfo(workId, edgeId, l, r);
    }

    // 分配资源
    for (int i = 0; i < changeEdgeInfoVec.size(); i ++ ) {
        int workId = changeEdgeInfoVec[i][0], edgeId = changeEdgeInfoVec[i][1], l = changeEdgeInfoVec[i][2], r = changeEdgeInfoVec[i][3];
        changeReInfo(workId, edgeId, l, r);
    }

    // 增加节点的变通道次数
    for (auto& [nodeId, workId] : addNodeAlternate) {
        p[nodeId].curUseP -= 1, p[nodeId].remainUseP += 1;
        nodeChangeHasWorkReal[nodeId][workId] = false;
        
    }
    // 减少节点的变通道次数
    for (auto& [nodeId, workId] : subNodeAlternate) {
        p[nodeId].curUseP += 1, p[nodeId].remainUseP -= 1;

        nodeChangeHasWorkReal[nodeId][workId] = true;
        nodeChangeHasWork[nodeId].insert(workId);
    }


    // //输出每个节点的变通道数量
    // for (int i = 1; i <= n; i ++ ) {
    //     cout << "node" << i << " " << p[i].remainUseP << endl;
    // }

}

void initFindPathDaTa(int workId) {
    memset(stt, 0, sizeof stt);
    memset(cnt, 0, sizeof cnt);
    memset(st, 0, sizeof st);
    memset(pre3, -1, sizeof pre3);
    memset(pre2, -1, sizeof pre2);
    memset(pre, -1, sizeof pre);
    for (int i = 1; i <= m; i ++ ) {
        memo[i] = {};
    }
    memset(memoSt, 0,sizeof memoSt);

    idx = 1;
    ttime = 1;
}

// 计算两个集合的Jaccard相似度系数
double jaccardSimilarity(const set<int>& set1, const set<int>& set2) {
    if (set1.size() == 0 || set2.size() == 0) return 0.0;

    set<int> intersection;
    set<int> unionSet;
    
    // 计算交集和并集
    set_intersection(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(intersection, intersection.begin()));
    set_union(set1.begin(), set1.end(), set2.begin(), set2.end(), inserter(unionSet, unionSet.begin()));

    // 计算Jaccard相似度系数
    return (1.0 * intersection.size()) / (1.0 * unionSet.size());
}

// 随机生成子集
std::set<int> generateSubset(const std::vector<int>& elements, int subsetSize, std::mt19937& gen) {
    std::vector<int> indices(elements.size());
    std::iota(indices.begin(), indices.end(), 0); // 生成索引

    // Fisher-Yates 洗牌算法随机选择子集
    for (int i = 0; i < subsetSize; ++i) {
        std::uniform_int_distribution<int> dist(i, elements.size() - 1);
        int j = dist(gen);
        std::swap(indices[i], indices[j]);
    }

    std::set<int> subset;
    for (int i = 0; i < subsetSize; ++i) {
        subset.insert(elements[indices[i]]);
    }

    return subset;
}
// 启发式搜索函数
std::vector<std::set<int>> heuristicSearch(const std::vector<std::set<int>>& sets, const std::set<int>& originalSet, int subsetSize, int numSubsets) {
    std::vector<std::set<int>> subsets;
    std::vector<int> elements(originalSet.begin(), originalSet.end());
    std::random_device rd;
    std::mt19937 gen(rd());
    int acc = 0;
    while (subsets.size() < numSubsets && acc <= 100000) {

        std::set<int> subset = generateSubset(elements, subsetSize, gen);
        bool isUnique = true;
        // 与原来的集合
        for (const auto& existingSubset : sets) {
            if (jaccardSimilarity(subset, existingSubset) > 0.5) {
                isUnique = false;
                break;
            }
        }
        // 与自己比较
        for (const auto& existingSubset : subsets) {
            if (jaccardSimilarity(subset, existingSubset) > 0.5) {
                isUnique = false;
                break;
            }
        }

        if (isUnique) {
            subsets.push_back(subset);
        }
        acc ++;
    }

    return subsets;
}

void dataProcessing() {
    set<int> originalSet; // 原始集合
    vector<PDI> sortG;

    for (int i = 1; i <= m; i ++ ) {
        auto [a, b] = edgeToNode[i];
        sortG.push_back({min(p[a].remainUseP, p[b].remainUseP), i});
    }
    sort(sortG.begin(), sortG.end());

    for (auto [_, id] : sortG) {
        originalSet.insert(id);
        if ((int)originalSet.size() >= 1.0 * m / 1.5){
            break;
        }
    }

    int m1 = originalSet.size();
    int subsetSize = min((int)(1.0 * m1 / 3), 60);
    int numSubsets = 30;
    vector<set<int>> subsets;
    while (subsetSize >= 1 && (int)subsets.size() < 30) {
        vector<set<int>> ssubsets = heuristicSearch(subsets, originalSet, subsetSize, numSubsets - subsets.size());
        subsetSize --;
        subsets.insert(subsets.end(), ssubsets.begin(), ssubsets.end());
    }

    // for (const auto& subset : subsets) {
    //     std::cout << "子集: ";
    //     for (const auto& element : subset) {
    //         std::cout << element << " ";
    //     }
    //     std::cout << std::endl;
    // }
    while (subsets.size() > 30) subsets.pop_back();
    saveBadEdge = subsets;

    T1 = saveBadEdge.size();
    // // 输出所有子集
    // for (int i = 0; i < subsets.size(); ++i) {
    //     std::cout << "Subset " << i + 1 << ": ";
    //     for (int num : subsets[i]) {
    //         std::cout << num << " ";
    //     }
    //     std::cout << std::endl;
    // }

}

void constructTest1() {

    cout << (int)saveBadEdge.size() << endl;
    fflush(stdout);

    for (int i = 0; i < saveBadEdge.size(); i ++ ) {
        // 输出每个场景的损坏的边的数量
        cout << saveBadEdge[i].size() << endl;
        fflush(stdout);
        for (int edgeId : saveBadEdge[i]) {
            cout << edgeId << " ";
        }
        cout << endl;
        fflush(stdout);
    }
}

void constructionPre3(int index, int workId) {
    int temp = index - 1;
    while (temp != 0) {
        int start = pre2[temp].x.x, end = pre2[temp].x.y, edgeId = pre2[temp].y.x, channelL = pre2[temp].y.y;
        // cout << start << "    " << end << " ";
        pre3[end] = {{start, edgeId}, {channelL, channelL + business[workId].width - 1}};
        temp = pre[temp];
    }
    // cout << endl;

}




