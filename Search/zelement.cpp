/*
zelement.cpp
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
#include "zelement.h"
#include <sys/time.h>
#include <immintrin.h>
#include <cstring>
#include <string>
#include "../utils.h"
#include <stdio.h>
#include <cmath>
#include <iostream>
#include "../Algorithm/zoxc.h"
extern "C"{
#include "../Jerasure-1.2A/jerasure.h"
#include "../Jerasure-1.2A/cauchy.h"
}
#include <cassert>

//改一下
//int ZElement::cpy_weight = 1532;
//int ZElement::xor_weight = 2999;

//测试了10次，每次10个循环，求平均值
int ZElement::cpy_weight = 1737;
int ZElement::xor_weight = 2104;

int ZElement::K = 0;
int ZElement::M = 0;
int ZElement::W = 0;
int ZElement::cost_func = 0;
int ZElement::strategy = 0;
bool ZElement::isInited = false;

FILE * ZElement::fp = NULL;

int ZElement::_chunk_size = 64;//MB
int ZElement::_packet_size = 64;//MB

int ZElement::_network_bandwidth = 1;//Gb/s
int ZElement::_disk_bandwidth = 100;//Mb/s

int ZElement::_chunk_num = 1000;//每个node的chunk数量
double ZElement::ec_transmit_pkt_time = 0;
double ZElement::ec_transmit_pkt_time2 = 0;
double ZElement::dir_transit_pkt_time = 0;
//_disk_bandwidth上下行相同？读取和写入时间相同？
double ZElement::rd_pkt_time = 0;
double ZElement::wr_pkt_time = 0;
double ZElement::repair_chunk_time = 0;
double ZElement::migrate_chunk_time = 0;
int    ZElement::migrateTotalCnt = 0;
double ZElement::totalERTime = 0;
double ZElement::totalSeanTime = 0;
double ZElement::totalERTempTime = 0;
/* 如果sch[idx][2] == 1|3|5，xor++ */
long long schedule_weight_3(vector<int*> &sch)
{
    long long ret = 0;
    int i;
    int size = sch.size();
    for(i = 0;i<size;i++)
    {
        switch(ZElement::cost_func)
        {
        case ZElement::COST_XOR:
            if(sch[i][2] == 1 || sch[i][2] == 3 || sch[i][2] == 5) ret++; // consider only xor
            break;
        case ZElement::COST_SUM:
            ret ++;
            break;
        case ZElement::COST_WEIGHTED:
            if(sch[i][2] == 1 || sch[i][2] == 3 || sch[i][2] == 5) // add weight to xor
                ret += ZElement::xor_weight;
            else
                ret += ZElement::cpy_weight;
            break;
        }
    }
    return ret;
}
/*如果sch[idx][4] == 1 ，xor++；*/
long long schedule_len_5(int** sch)
{
    long long ret = 0;
    int idx = 0;
    int cpCount = 0;
    while(sch[idx][0] != -1)
    {
        switch(ZElement::cost_func)
        {
        case ZElement::COST_XOR:
            if(sch[idx][4] == 1) ret++; // consider only xor
            break;
        case ZElement::COST_SUM:
            ret ++;
            break;
        case ZElement::COST_WEIGHTED:
            if(sch[idx][4] == 1) // add weight to xor and memcpy
                ret += ZElement::xor_weight;
            else{
                cpCount++;
                ret += ZElement::cpy_weight;
            }

            break;
        }
        idx ++;
    }
    printf("the schedule cpCount = %d\n",cpCount);
    return ret;
}

ZElement::ZElement(int *p)
{
    assert(isInited);
    array = vector<int>(p,p+K+M); //vector初值为p+K+M


}

ZElement::ZElement(vector<int> p)
{
    assert(isInited);
    array = p;
}
// calculate the number of chunks to be migrated based on the mathematical model
int ZElement::calMigrateChunkNum(int num_repair_chunk){

    int num_repair_chunk_node = 1;//需要修复节点数



//    cout << "num_repair_chunk_node = " << num_repair_chunk_node << endl;
    // for case of pipelined repair
//    double repair_chunk_time;
//    double migrate_chunk_time;
    int packet_num;

    packet_num = _chunk_size/_packet_size;
    // if transmission time > read/write time
    //_packet_size(MB)，_network_bandwidth(Gb/s),time(s)
    ec_transmit_pkt_time = _packet_size*K*num_repair_chunk_node*8*1.0/(1024*_network_bandwidth);
    ec_transmit_pkt_time2 = _packet_size*(K-1)*num_repair_chunk_node*8*1.0/(1024*_network_bandwidth);//在ER->H少传一份数据
    dir_transit_pkt_time = _packet_size*8*1.0/(1024*_network_bandwidth);
    //_disk_bandwidth上下行相同？读取和写入时间相同？
    rd_pkt_time = _packet_size*1.0/_disk_bandwidth;
    wr_pkt_time = _packet_size*1.0/_disk_bandwidth;

//    printf("_chunk_size:%d, _packet_size:%d, K:%d, _network_bandwidth:%d, _disk_bandwidth:%d\n",_chunk_size,_packet_size,
//           K,_network_bandwidth,_disk_bandwidth);

//    printf("ec_transmit_pkt_time = %f, dir_transit_pkt_time = %f\n",ec_transmit_pkt_time,dir_transit_pkt_time);


    // determine the repair time
    if(ec_transmit_pkt_time > wr_pkt_time)
        repair_chunk_time = ec_transmit_pkt_time*packet_num + wr_pkt_time*num_repair_chunk_node + rd_pkt_time;
    else
        repair_chunk_time = wr_pkt_time*packet_num*num_repair_chunk_node + ec_transmit_pkt_time + rd_pkt_time;

    // determine the migration time
    if(dir_transit_pkt_time > wr_pkt_time)
        migrate_chunk_time = dir_transit_pkt_time*packet_num + wr_pkt_time + rd_pkt_time;
    else
        migrate_chunk_time = wr_pkt_time*packet_num + dir_transit_pkt_time + rd_pkt_time;

//    cout << "repair_chunk_time = " << repair_chunk_time << endl;
//    cout << "migrate_chunk_time = " << migrate_chunk_time << endl;

    int num_least_migrate_chunk = repair_chunk_time/migrate_chunk_time;
    double normal_aver_time = repair_chunk_time/(num_repair_chunk+num_least_migrate_chunk);
    double pref_mgrt_aver_time = (num_least_migrate_chunk+1)*migrate_chunk_time/(num_repair_chunk+num_least_migrate_chunk+1);

    if(normal_aver_time < pref_mgrt_aver_time) return num_least_migrate_chunk;
    else return num_least_migrate_chunk+1;

}

// calculate the number of chunks to be migrated based on the mathematical model
// 第二阶段，将H作为新数据盘，N作为新的H（parity node）
int ZElement::calMigrateChunkNum2(int num_repair_chunk){

    int num_repair_chunk_node = 1;//需要修复节点数

    int packet_num;

    packet_num = _chunk_size/_packet_size;
    // if transmission time > read/write time
    //_packet_size(MB)，_network_bandwidth(Gb/s),time(s)
//    ec_transmit_pkt_time = _packet_size*K*num_repair_chunk_node*8*1.0/(1024*_network_bandwidth);
//    dir_transit_pkt_time = _packet_size*8*1.0/(1024*_network_bandwidth);
//    //_disk_bandwidth上下行相同？读取和写入时间相同？
//    rd_pkt_time = _packet_size*1.0/_disk_bandwidth;
//    wr_pkt_time = _packet_size*1.0/_disk_bandwidth;

//    printf("_chunk_size:%d, _packet_size:%d, K:%d, _network_bandwidth:%d, _disk_bandwidth:%d\n",_chunk_size,_packet_size,
//           K,_network_bandwidth,_disk_bandwidth);

//    printf("ec_transmit_pkt_time = %f, dir_transit_pkt_time = %f\n",ec_transmit_pkt_time,dir_transit_pkt_time);


    // determine the repair time
//    if(ec_transmit_pkt_time > wr_pkt_time)
//        repair_chunk_time = ec_transmit_pkt_time*packet_num + wr_pkt_time*num_repair_chunk_node + rd_pkt_time;
//    else
//        repair_chunk_time = wr_pkt_time*packet_num*num_repair_chunk_node + ec_transmit_pkt_time + rd_pkt_time;

//    // determine the migration time
//    if(dir_transit_pkt_time > wr_pkt_time)
//        migrate_chunk_time = dir_transit_pkt_time*packet_num + wr_pkt_time + rd_pkt_time;
//    else
//        migrate_chunk_time = wr_pkt_time*packet_num + dir_transit_pkt_time + rd_pkt_time;

//    cout << "repair_chunk_time = " << repair_chunk_time << endl;
//    cout << "migrate_chunk_time = " << migrate_chunk_time << endl;



    int num_least_migrate_chunk = repair_chunk_time/migrate_chunk_time;
    double normal_aver_time = repair_chunk_time/(num_repair_chunk+num_least_migrate_chunk);
    double pref_mgrt_aver_time = (num_least_migrate_chunk+1)*migrate_chunk_time/(num_repair_chunk+num_least_migrate_chunk+1);

    if(normal_aver_time < pref_mgrt_aver_time) return num_least_migrate_chunk;
    else return num_least_migrate_chunk+1;

}
int ZElement::doRepair(queue<int> & repairMatrix, queue<int> & migrateMatrix, double limitTime){

    double curTime = 0;
    int repairCnt = 0;
    int migrateCnt = 0;
    int round = 0;
    while(true){

        if(curTime > limitTime) break;
        if(repairMatrix.size() == 0 && migrateMatrix.size() == 0) break;

        printf("round %d, time:%lf \n",round,curTime);
        curTime += repair_chunk_time;

        migrateCnt = calMigrateChunkNum(1);


        if(repairMatrix.size() > 0){
            printf("修复块：%d ",repairMatrix.front());
            repairMatrix.pop();
        }
        if(migrateMatrix.size() > 0)
            printf("迁移块：");
        for(int i=0;i<migrateCnt && migrateMatrix.size() > 0;i++) {
            printf("%d ",migrateMatrix.front());
            migrateMatrix.pop();
        }
        printf("\n");

        repairCnt+=migrateCnt + 1;
        migrateTotalCnt+=migrateCnt;

        round++;
    }

    return repairCnt;
}
//执行第一阶段：完全并发
int ZElement::doRepairBeforeConcurrent(int U, int P, double limitTime){

    double migrateTime = 0;//迁移U-P所需时间
    double restTime = 0;//迁移完剩余时间
    double repairTime = 0;//实际的修复所花时间
    int packetNum = _chunk_size / _packet_size;
    int migrateCnt = 0, repairCnt = 0, migrRealCnt = 0;
    migrRealCnt = U - P;

    //迁移
    if(migrRealCnt > 0){

        if(rd_pkt_time > dir_transit_pkt_time){
            migrateCnt = ( limitTime - dir_transit_pkt_time - wr_pkt_time) / (wr_pkt_time * packetNum) ;
            migrateTime = (U - P)*rd_pkt_time*packetNum + dir_transit_pkt_time + wr_pkt_time;
        }else{
            migrateCnt = ( limitTime - rd_pkt_time - wr_pkt_time) / (dir_transit_pkt_time * packetNum) ;
            migrateTime = (U - P)*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        //limitTime < migrateTime
        if( migrateCnt < ( U - P ) ){
            migrRealCnt = migrateCnt;
            printf("经过%lfs，已迁移：%d-%d,共%d块\n",limitTime,P,P+migrRealCnt-1,migrRealCnt);
        }
        else
            printf("经过%lfs，已迁移：%d-%d,共%d块\n",migrateTime,P,P+migrRealCnt-1,migrRealCnt);
    }
    //修复
    if(wr_pkt_time > ec_transmit_pkt_time2){
        repairCnt = ( limitTime - ec_transmit_pkt_time2 - wr_pkt_time) / (wr_pkt_time * packetNum)  ;
        repairTime = P*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
    }else{
        repairCnt = ( limitTime - rd_pkt_time - wr_pkt_time) / (ec_transmit_pkt_time2 * packetNum)  ;
        repairTime = P*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
    }
    repairTime = repairTime < limitTime ? repairTime : limitTime;
    repairCnt = repairCnt < P ? repairCnt : P;
    if(repairCnt > 0){
        printf("经过%lfs，已修复：%d-%d,共%d块\n",repairTime,0,repairCnt-1,repairCnt);
    }



    totalSeanTime += repairTime;

    printf("剩余: %d块需要修复！\n",U - repairCnt - migrRealCnt);

    return repairCnt + migrRealCnt;
}
int ZElement::doRepairAfterConcurrent(int U, int P,int schedule){

    double migrateTime = 0;//迁移U-P所需时间
    double restTime = 0;//迁移完剩余时间
    double er2hTime = 0;//实际的修复所花时间
    double er2nTime = 0;
    double T1 = 0,T2 = 0;
    int packetNum = _chunk_size / _packet_size;
    int migrateCnt = 0, repairCnt = 0, migrRealCnt = 0;
    double migrateSinglePacketTime = 0, repairSinglePacketTime =0;
//    int U = 0,P = 0;
//    P = repairMatrix.size();
//    U = P + migrateMatrix.size();

    switch (schedule) {
    //H作为数据盘
    case 0:
        //迁移：H->N
        if( (U - P) > 0 ){
            if(rd_pkt_time > dir_transit_pkt_time){
                migrateTime = (U - P)*rd_pkt_time*packetNum + dir_transit_pkt_time + wr_pkt_time;
            }else{
                migrateTime = (U - P)*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
            }
            //修复：ER->H，H已有一份数据，用ec_transmit_pkt_time2
            if(rd_pkt_time > ec_transmit_pkt_time2){
                er2hTime = (U - P)*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
            }else{
                er2hTime = (U - P)*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
            }
            T1 = max(migrateTime,er2hTime);
            printf("并发完成修复：ER->H, 迁移：H->N %d块 工作，耗时：%lfs\n",U-P,T1);
        }
        //修复：ER->N
        if(rd_pkt_time > ec_transmit_pkt_time){
            er2nTime = P*rd_pkt_time*packetNum + ec_transmit_pkt_time + wr_pkt_time;
        }else{
            er2nTime = P*ec_transmit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T2 = er2nTime;


        totalSeanTime += max(T1,T2);

        migrateSinglePacketTime = rd_pkt_time+dir_transit_pkt_time+wr_pkt_time;
        repairSinglePacketTime = rd_pkt_time+ec_transmit_pkt_time2+wr_pkt_time;
        printf("H->N迁移单包时间：%lfs,读：%lfs,传：%lfs,写：%lfs\n",migrateSinglePacketTime,rd_pkt_time,dir_transit_pkt_time,wr_pkt_time);
        printf("H->N迁移单包时间：%lfs,读：%lfs,传：%lfs,写：%lfs\n",repairSinglePacketTime,rd_pkt_time,ec_transmit_pkt_time2,wr_pkt_time);

        printf("并发完成修复：ER->N %d块 工作，耗时：%lfs\n",P,T2);
        printf("将H作为数据盘，N作为parity，总耗时：%lfs\n",T1+T2);
        break;
    //H作为parity
    case 1:
        //迁移：H->N，将H中的data迁移给N
        if( P > 0 ){
            if(wr_pkt_time > dir_transit_pkt_time){
                migrateTime = P*wr_pkt_time*packetNum + dir_transit_pkt_time + rd_pkt_time;
            }else{
                migrateTime = P*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
            }

            //修复：ER->H，H已有一份数据，用ec_transmit_pkt_time2
            if(wr_pkt_time > ec_transmit_pkt_time2){
                er2hTime = P*wr_pkt_time*packetNum + ec_transmit_pkt_time2 + rd_pkt_time;
            }else{
                er2hTime = P*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
            }
            T1 = max(migrateTime,er2hTime);
            printf("并发完成修复：ER->H, 迁移：H->N %d块 工作，耗时：%lfs\n",U-P,T1);
        }
        //修复：ER->N
        if(wr_pkt_time > ec_transmit_pkt_time){
            er2nTime = (U-P)*wr_pkt_time*packetNum + ec_transmit_pkt_time + rd_pkt_time;
        }else{
            er2nTime = (U-P)*ec_transmit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T2 = er2nTime;

        totalSeanTime += max(T1,T2);

        printf("并发完成修复：ER->N %d块 工作，耗时：%lfs\n",P,T2);
        printf("将H作为数据盘，N作为parity，总耗时：%lfs\n",T1+T2);
        break;
    default:
        break;
    }


    return 0;
}
//执行第一阶段
int ZElement::doRepairBeforeNewDisk(int U, int P, double limitTime){

    double migrateTime = 0;//迁移U-P所需时间
    double restTime = 0;//迁移完剩余时间
    double repairTime = 0;//实际的修复所花时间
    int packetNum = _chunk_size / _packet_size;
    int migrateCnt = 0, repairCnt = 0, migrRealCnt = 0;
    migrRealCnt = U - P;

    //迁移
    if(rd_pkt_time > dir_transit_pkt_time){
        migrateCnt = ( limitTime - dir_transit_pkt_time - wr_pkt_time) / (rd_pkt_time * packetNum) ;
        migrateTime = (U - P)*rd_pkt_time*packetNum + dir_transit_pkt_time + wr_pkt_time;
    }else{
        migrateCnt = ( limitTime - rd_pkt_time - wr_pkt_time) / (dir_transit_pkt_time * packetNum) ;
        migrateTime = (U - P)*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
    }
    //limitTime < migrateTime
    if( migrateCnt < ( U - P ) ){
        migrRealCnt = migrateCnt;
        printf("经过%lfs，已迁移：%d-%d,共%d块\n",limitTime,P,P+migrRealCnt-1,migrRealCnt);
    }
    else{
        printf("经过%lfs，已迁移：%d-%d,共%d块\n",migrateTime,P,P+migrRealCnt-1,migrRealCnt);
        //修复
        restTime = limitTime - migrateTime;
        if(rd_pkt_time > ec_transmit_pkt_time2){
            repairCnt = ( restTime - ec_transmit_pkt_time2 - wr_pkt_time) / (rd_pkt_time * packetNum)  ;
            repairTime = P*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
        }else{
            repairCnt = ( restTime - rd_pkt_time - wr_pkt_time) / (ec_transmit_pkt_time2 * packetNum)  ;
            repairTime = P*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
        }
        repairTime = repairTime < restTime ? repairTime : restTime;
        repairCnt = repairCnt < P ? repairCnt : P;
        if(repairCnt > 0){
            printf("经过%lfs，已修复：%d-%d,共%d块\n",repairTime,0,repairCnt-1,repairCnt);
        }

    }

    totalSeanTime += limitTime;

    printf("剩余: %d块需要修复！\n",U - repairCnt - migrRealCnt);

    return repairCnt + migrRealCnt;
}
//执行第一阶段
int ZElement::doRepairBeforeOnlyERTemp(int U){

    double repairTime = 0;//实际的修复所花时间
    int packetNum = _chunk_size / _packet_size;

    //修复
    if(rd_pkt_time > ec_transmit_pkt_time){
        repairTime = U*rd_pkt_time*packetNum + ec_transmit_pkt_time ;
    }else{
        repairTime = U*ec_transmit_pkt_time*packetNum + rd_pkt_time ;
    }

    totalERTempTime += repairTime;

    printf("经过%lfs，已修复：%d-%d,共%d块\n",repairTime,0,U-1,U);
    return U;
}

//执行第一阶段
int ZElement::doRepairBeforeNewDiskOnlyER(int U){

    double repairTime = 0;//实际的修复所花时间
    int packetNum = _chunk_size / _packet_size;

    //修复
    if(rd_pkt_time > ec_transmit_pkt_time2){
        repairTime = U*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
    }else{
        repairTime = U*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
    }

    totalERTime += repairTime;

    printf("经过%lfs，已修复：%d-%d,共%d块\n",repairTime,0,U-1,U);
    return U;
}
//执行第二阶段
double ZElement::doRepairAfterNewDiskOnlyER(int U){

    double repairTime = 0;//实际的修复所花时间
    int packetNum = _chunk_size / _packet_size;

    //修复
    if(rd_pkt_time > ec_transmit_pkt_time){
        repairTime = U*rd_pkt_time*packetNum + ec_transmit_pkt_time + wr_pkt_time;
    }else{
        repairTime = U*ec_transmit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
    }


    printf("经过%lfs，已修复：%d-%d,共%d块\n",repairTime,0,U-1,U);
    return repairTime;
}
int ZElement::doRepairAfterNewDisk(int U, int P,int schedule){

    double migrateTime = 0;//迁移U-P所需时间
    double restTime = 0;//迁移完剩余时间
    double er2hTime = 0;//实际的修复所花时间
    double er2nTime = 0;
    double T1 = 0,T2 = 0;
    int packetNum = _chunk_size / _packet_size;
    int migrateCnt = 0, repairCnt = 0, migrRealCnt = 0;
    double migrateSinglePacketTime = 0, repairSinglePacketTime =0;
//    int U = 0,P = 0;
//    P = repairMatrix.size();
//    U = P + migrateMatrix.size();

    switch (schedule) {
    //H作为数据盘
    case 0:
        //迁移：H->N
        if(rd_pkt_time > dir_transit_pkt_time){
            migrateTime = (U - P)*rd_pkt_time*packetNum + dir_transit_pkt_time + wr_pkt_time;
        }else{
            migrateTime = (U - P)*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        //修复：ER->H，H已有一份数据，用ec_transmit_pkt_time2
        if(rd_pkt_time > ec_transmit_pkt_time2){
            er2hTime = (U - P)*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
        }else{
            er2hTime = (U - P)*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T1 = migrateTime+er2hTime;
        printf("顺序完成修复：ER->H, 迁移：H->N %d块 工作，耗时：%lfs\n",U-P,T1);

        //修复：ER->N
        if(rd_pkt_time > ec_transmit_pkt_time){
            er2nTime = P*rd_pkt_time*packetNum + ec_transmit_pkt_time + wr_pkt_time;
        }else{
            er2nTime = P*ec_transmit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T2 = er2nTime;

        totalSeanTime += T1 + T2;

        migrateSinglePacketTime = rd_pkt_time+dir_transit_pkt_time+wr_pkt_time;
        repairSinglePacketTime = rd_pkt_time+ec_transmit_pkt_time2+wr_pkt_time;
        printf("H->N迁移单包时间：%lfs,读：%lfs,传：%lfs,写：%lfs\n",migrateSinglePacketTime,rd_pkt_time,dir_transit_pkt_time,wr_pkt_time);
        printf("H->N迁移单包时间：%lfs,读：%lfs,传：%lfs,写：%lfs\n",repairSinglePacketTime,rd_pkt_time,ec_transmit_pkt_time2,wr_pkt_time);

        printf("完成修复：ER->N %d块 工作，耗时：%lfs\n",P,T2);
        printf("将H作为数据盘，N作为parity，总耗时：%lfs\n",T1+T2);
        break;
    //H作为parity
    case 1:
        //迁移：H->N，将H中的data迁移给N
        if(rd_pkt_time > dir_transit_pkt_time){
            migrateTime = P*rd_pkt_time*packetNum + dir_transit_pkt_time + wr_pkt_time;
        }else{
            migrateTime = P*dir_transit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        //修复：ER->H，H已有一份数据，用ec_transmit_pkt_time2
        if(rd_pkt_time > ec_transmit_pkt_time2){
            er2hTime = P*rd_pkt_time*packetNum + ec_transmit_pkt_time2 + wr_pkt_time;
        }else{
            er2hTime = P*ec_transmit_pkt_time2*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T1 = migrateTime+er2hTime;
        printf("同时完成修复：ER->H, 迁移：H->N %d块 工作，耗时：%lfs\n",U-P,T1);

        //修复：ER->N
        if(rd_pkt_time > ec_transmit_pkt_time){
            er2nTime = (U-P)*rd_pkt_time*packetNum + ec_transmit_pkt_time + wr_pkt_time;
        }else{
            er2nTime = (U-P)*ec_transmit_pkt_time*packetNum + rd_pkt_time + wr_pkt_time;
        }
        T2 = er2nTime;

        totalSeanTime += T1 + T2;

        printf("完成修复：ER->N %d块 工作，耗时：%lfs\n",P,T2);
        printf("将H作为数据盘，N作为parity，总耗时：%lfs\n",T1+T2);
        break;
    default:
        break;
    }


    return 0;
}


void ZElement::init(int tK, int tM, int tW, int tcost, int tstrategy)
{
    K = tK;
    M = tM;
    W = tW;
    cost_func = tcost;
    strategy = tstrategy;
    isInited = true;
}

long long ZElement::value()
{
    long long xor0_natual_dumb;
    long long my_xor0_natual_dumb;
    long long xor1_natual_smart;
    long long xor2_natual_grouping_unweighted;
    long long xor3_natual_grouping_weighted;
    long long xor4_normal_dumb;
    long long xor5_normal_smart;
    long long xor6_normal_grouping_unweighted;
    long long xor7_normal_grouping_weighted;

    ZOXC *xc;
    int *matrix;
    int *bitmatrix;
    int **schedule;
    int **myschedule;
    //X:array[K----K+m];Y:array[0----K]
    matrix = cauchy_xy_coding_matrix(K,M,W,array.data()+K,array.data());


    switch(strategy)
    {
    case 0:

        // natual dumb
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        jerasure_print_bitmatrix(bitmatrix,K*W,M*W,W);
        //myschedule = arr_trace(K,M,W,bitmatrix);
        schedule = jerasure_dumb_bitmatrix_to_schedule(K,M,W,bitmatrix);
        xor0_natual_dumb = schedule_len_5(schedule); //计算XOR的cost
        my_xor0_natual_dumb = schedule_len_5(schedule); //计算XOR的cost

        fprintf(fp,"组(%d,%d,%d)结果不同：原cost：%d,现在cost：%d\n",
                                                           K,M,W,xor0_natual_dumb,my_xor0_natual_dumb);
        free2dSchedule5(schedule); 
        free(bitmatrix);
        free(matrix);
        v = xor0_natual_dumb;
        fclose(fp);
//        free(fp);
        return xor0_natual_dumb;
    case 1:
        // natual smart
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        schedule = jerasure_smart_bitmatrix_to_schedule(K,M,W,bitmatrix);
        xor1_natual_smart = schedule_len_5(schedule);
        free2dSchedule5(schedule);
        free(bitmatrix);
        free(matrix);
        v = xor1_natual_smart;
        return xor1_natual_smart;
    case 2:
        // natual grouping unweighed
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        xc = new ZOXC(K*W, M*W);
        xc->grouping_1s(bitmatrix, false);
        xor2_natual_grouping_unweighted = schedule_weight_3(xc->schedule);
        free(bitmatrix);
        free(matrix);
        delete xc;
        v = xor2_natual_grouping_unweighted;
        return xor2_natual_grouping_unweighted;
    case 3:
        // natual grouping weightd
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        xc = new ZOXC(K*W, M*W);
        xc->grouping_1s(bitmatrix, true);
        xor3_natual_grouping_weighted = schedule_weight_3(xc->schedule);
        free(bitmatrix);
        free(matrix);
        delete xc;
        v = xor3_natual_grouping_weighted;
        return xor3_natual_grouping_weighted;
    case 4:
        // normal dumb
        cauchy_improve_coding_matrix(K,M,W,matrix);
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        schedule = jerasure_dumb_bitmatrix_to_schedule(K,M,W,bitmatrix);
        xor4_normal_dumb = schedule_len_5(schedule);
        free2dSchedule5(schedule);
        free(bitmatrix);
        free(matrix);
        v = xor4_normal_dumb;
        return xor4_normal_dumb;
    case 5:
        // normal smart
        cauchy_improve_coding_matrix(K,M,W,matrix);
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        schedule = jerasure_smart_bitmatrix_to_schedule(K,M,W,bitmatrix);
        xor5_normal_smart = schedule_len_5(schedule);
        free2dSchedule5(schedule);
        free(bitmatrix);
        free(matrix);
        v = xor5_normal_smart;
        return xor5_normal_smart;
    case 6:
        // normal unweighted
        cauchy_improve_coding_matrix(K,M,W,matrix);
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        xc = new ZOXC(K*W, M*W);
        xc->grouping_1s(bitmatrix, false);
        xor6_normal_grouping_unweighted = schedule_weight_3(xc->schedule);
        free(bitmatrix);
        free(matrix);
        delete xc;
        v = xor6_normal_grouping_unweighted;
        return xor6_normal_grouping_unweighted;
    case 7:
        // normal weighted
        cauchy_improve_coding_matrix(K,M,W,matrix);
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
//        jerasure_print_matrix(matrix,M,K,W);
//        jerasure_print_bitmatrix(bitmatrix,M*W,K*W,W);
        xc = new ZOXC(K*W, M*W);
        xc->grouping_1s(bitmatrix, true);
        xor7_normal_grouping_weighted = schedule_weight_3(xc->schedule);
        free(bitmatrix);
        free(matrix);
        delete xc;
        v = xor7_normal_grouping_weighted;
        return xor7_normal_grouping_weighted;

    // Sean strategy
    case 8:

        double p = Random(1,10)*1.0/10;//产生10种概率，分别为0.1，0.2,...,1
        p = 1 - 0.3;
        int U = _chunk_num;//需要修复的块的总数量
        int P = U * p ;//第一阶段坏的块的数量
//        int migrBestCnt = calMigrateChunkNum(U);
//        int migrRealCnt = U - P;
        int mt = Random(1,1000)*1.0;//随机产生失效时间
        mt = 500;
        int i = 0,repairCnt = 0,restCnt = U,restRepairCnt = 0,migrateCnt;
//        double migrateTime = 0,restTime = 0,repairTime = 0;
//        int packetSize = _chunk_size / _packet_size;

        double resERTime = 0, resSeanTime = 0;
        queue<int> migrateMatrix,repairMatrix;
        for(; i < P ;i++){
            repairMatrix.push(i);
        }
        for(;i < U ;i++){
            migrateMatrix.push(i);
        }
        //初始化数据
        ec_transmit_pkt_time = _packet_size*K*8*1.0/(1024*_network_bandwidth);
        ec_transmit_pkt_time2 = _packet_size*(K-1)*8*1.0/(1024*_network_bandwidth);
        dir_transit_pkt_time = _packet_size*8*1.0/(1024*_network_bandwidth);
        rd_pkt_time = _packet_size*1.0/_disk_bandwidth;
        wr_pkt_time = _packet_size*1.0/_disk_bandwidth;
        int packetNum = _chunk_size / _packet_size;
        if(rd_pkt_time > ec_transmit_pkt_time){
            resERTime = packetNum*rd_pkt_time + ec_transmit_pkt_time ;
        }else{
            resERTime = packetNum*ec_transmit_pkt_time + rd_pkt_time ;
        }

        if(rd_pkt_time > ec_transmit_pkt_time2){
            resSeanTime = packetNum*rd_pkt_time + ec_transmit_pkt_time2 ;
        }else{
            resSeanTime = packetNum*ec_transmit_pkt_time2 + rd_pkt_time ;
        }






        printf("************仿真开始*************\n");

        printf("*******纯ER算法（保证持久化）:开始*******\n");
        printf("第一阶段开始:\n");
        doRepairBeforeNewDiskOnlyER(U);
        printf("第二阶段开始:\n");
        totalERTime += doRepairAfterNewDiskOnlyER(U);
        printf("纯ER算法总耗时：%lf\n",totalERTime);
        printf("*******纯ER算法（保证持久化）:结束*******\n");


        printf("*******纯ER算法（非持久化）:开始*******\n");
        printf("第一阶段开始:\n");
        doRepairBeforeOnlyERTemp(P);
        printf("第二阶段开始:\n");
        totalERTempTime += doRepairAfterNewDiskOnlyER(U);
        printf("纯ER算法总耗时：%lf\n",totalERTempTime);
        printf("*******纯ER算法（非持久化）:结束*******\n");


        printf("*******并发Sean算法:开始*******\n");
        printf("第一阶段开始:\n");
        printf("失效数据量:%d, 在%ds之后，整个盘将失效...\n",P,mt);
        printf("ER恢复单包时间：%lf, 迁移单包时间：%lf\n",ec_transmit_pkt_time,dir_transit_pkt_time);

        restCnt = U - doRepairBeforeConcurrent(U,P,mt);


        printf("第二阶段开始:\n");
        if(restCnt > 0.5*U ){
            doRepairAfterConcurrent(U,U - restCnt, 1);
        }else{
            doRepairAfterConcurrent(U,U - restCnt, 0);
        }



        printf("Sean算法总耗时：%lf\n",totalSeanTime);
        printf("*******Sean算法:结束*******\n");

        printf("Sean算法：纯ER算法 总耗时对比 %.2f\n",totalSeanTime/totalERTime);
        printf("Sean算法：纯ER算法 数据响应时间对比  %.2f : %.2f\n",resSeanTime,resERTime);


        printf("*****************仿真结束*****************\n");
    }

    return 0;
}

// test the weight of memcpy and xor using given size of data, 1GB by default
// test 10 times, in each test run memcpy and xor *loops* times
void ZElement::test_cost_weight(int size, int loops)
{
    printf("Test_cost_weight: size %d, loops %d\n", size, loops);
    cpy_weight = 0;
    xor_weight = 0;
    for(int i = 0;i<10;i++)
    {
        char *dat, *dat1, *dat2;

        posix_memalign((void**)&dat, 32, size);
        posix_memalign((void**)&dat1, 32, size);
        posix_memalign((void**)&dat2, 32, size);

        /*
            * struct  timeval{
            *  long  tv_sec;/*秒*
            *  long  tv_usec;/*微妙*
            * }；
        */
        struct timeval t0,t1;

        gettimeofday(&t0,NULL);
        for(int j = 0;j<loops;j++)
            memcpy(dat1,dat,size);
        gettimeofday(&t1,NULL);

        cpy_weight += diff_us(t0,t1);

        gettimeofday(&t0,NULL);
        for(int j = 0;j<loops;j++)
            fast_xor(dat,dat1,dat1, size);
        gettimeofday(&t1,NULL);
        xor_weight += diff_us(t0,t1);

        free(dat);
        free(dat1);
        free(dat2);
    }
    printf("cpy_weight = %d, xor_weight = %d\n", cpy_weight, xor_weight);

    // scale to 0~9999 to avoid overflow
    while(cpy_weight > 10000)
    {
        cpy_weight /= 10;
        xor_weight /= 10;
    }
    if((fp = fopen("/Users/flyln/Documents/workplace/606/git/build-zerasure-Desktop_Qt_5_8_0_clang_64bit-Debug/result.txt","a+")) == NULL){
        printf("打开文件失败!\n");
        exit(-1);
    }

    fprintf(fp,"%d %d\n",cpy_weight,xor_weight);
    fclose(fp);
    printf("cpy_weight = %d, xor_weight = %d\n", cpy_weight, xor_weight);
}
