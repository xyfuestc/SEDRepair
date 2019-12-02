/*
zelement.h
Tianli Zhou

Fast Erasure Coding for Data Storage: A Comprehensive Study of the Acceleration Techniques

Revision 1.0
Mar 20, 2019

Tianli Zhou
Department of Electrical & Computer Engineering
Texas A&M University
College Station, TX, 77843
zhoutianli01@tamu.edu

Copyright (c) 2019, Tianli Zhou
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in
  the documentation and/or other materials provided with the
  distribution.

- Neither the name of the Texas A&M University nor the names of its
  contributors may be used to endorse or promote products derived
  from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ZELEMENT_H
#define ZELEMENT_H

#include <vector>
#include <queue>
#include <stdlib.h>

using namespace std;

class ZElement
{
public:
    enum cost_type{COST_XOR, COST_SUM, COST_WEIGHTED};
    ZElement(int *p=NULL);
    ZElement(vector<int> p);
    long long value();
    long long v;            // save the ret of value(), avoid re-compute

    static void test_cost_weight(int size = 100*1024*1024, int loops = 10);
    static int cpy_weight;
    static int xor_weight;

    static int _chunk_size;
    static int _packet_size;

    static int _network_bandwidth;
    static int _disk_bandwidth;
    static int _chunk_num;
    //单个packet的ec传输时间
    static double ec_transmit_pkt_time;
    static double ec_transmit_pkt_time2;
    static double dir_transit_pkt_time;
    //_disk_bandwidth上下行相同？读取和写入时间相同？
    static double rd_pkt_time;
    static double wr_pkt_time;
    static double repair_chunk_time;
    static double migrate_chunk_time;
    static int migrateTotalCnt;
    //总耗时对比
    static double totalERTime;
    static double totalERTempTime;
    static double totalSeanTime;

    // must be called before create any object
    static void init(int tK, int tM, int tW, int tcost, int tstrategy);

    static int K;
    static int M;
    static int W;
    static int cost_func; // 0: # of XOR, 1: # of ops, 2: weight sum
    static int strategy; // see paper, 0-7, strategy-(@/4,@%4)
    static bool isInited; // ensure init() has been called

    vector<int> array;

    //结果文件(result.txt)句柄
    static FILE *fp;

    int calMigrateChunkNum(int num_repair_chunk);
    int calMigrateChunkNum2(int num_repair_chunk);
    int doRepair(queue<int> & repairMatrix, queue<int> & migrateMatrix, double limitTime);
    int doRepairBeforeNewDisk(int U, int P, double limitTime);
    int doRepairBeforeNewDiskOnlyER(int U);
    double doRepairAfterNewDiskOnlyER(int U);
    int doRepairAfterNewDisk(int U, int P,int schedule);

    int doRepairBeforeConcurrent(int U, int P, double limitTime);
    int doRepairAfterConcurrent(int U, int P,int schedule);
    int doRepairBeforeOnlyERTemp(int U);
};

#endif // ZELEMENT_H
