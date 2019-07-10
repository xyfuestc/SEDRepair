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
#include "../utils.h"
#include <stdio.h>
#include "../Algorithm/zoxc.h"
extern "C"{
#include "../Jerasure-1.2A/jerasure.h"
#include "../Jerasure-1.2A/cauchy.h"
}
#include <cassert>

int ZElement::cpy_weight = 1532;
int ZElement::xor_weight = 2999;
int ZElement::K = 0;
int ZElement::M = 0;
int ZElement::W = 0;
int ZElement::cost_func = 0;
int ZElement::strategy = 0;
bool ZElement::isInited = false;
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
    FILE *fp = NULL;
    if((fp = fopen("/Users/flyln/Documents/workplace/606/git/build-zerasure-Desktop_Qt_5_8_0_clang_64bit-Debug/result.txt","a+")) == NULL){
        printf("打开文件失败!\n");
        exit(-1);
    }

    switch(strategy)
    {
    case 0:

        // natual dumb
        bitmatrix = jerasure_matrix_to_bitmatrix(K,M,W,matrix);
        jerasure_print_bitmatrix(bitmatrix,K*W,M*W,W);
        myschedule = arr_trace(K,M,W,bitmatrix);
        schedule = jerasure_dumb_bitmatrix_to_schedule(K,M,W,bitmatrix);
        xor0_natual_dumb = schedule_len_5(schedule); //计算XOR的cost
        my_xor0_natual_dumb = schedule_len_5(schedule); //计算XOR的cost

        if(xor0_natual_dumb != my_xor0_natual_dumb) fprintf(fp,"组(%d,%d,%d)结果不同：原cost：%d,现在cost：%d\n",
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
    printf("cpy_weight = %d, xor_weight = %d\n", cpy_weight, xor_weight);
}
