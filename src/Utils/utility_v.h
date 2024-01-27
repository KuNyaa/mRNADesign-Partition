#ifndef FASTCKY_UTILITY_V_H
#define FASTCKY_UTILITY_V_H

#include <string.h>
#include <cmath>
#include <assert.h>

#include "energy_parameter.h"
#include "intl11.h"
#include "intl21.h"
#include "intl22.h"

#define NUC_TO_PAIR(x,y) (x==1? (y==4?5:0) : (x==2? (y==3?1:0) : (x==3 ? (y==2?2:(y==4?3:0)) : (x==4 ? (y==3?4:(y==1?6:0)) : 0))))
#define PAIR_TO_LEFT_NUC(x) (x==1? 2:((x==2 || x==3)? 3:(x==5)? 1:4))
#define PAIR_TO_RIGHT_NUC(x) (x==2? 2:((x==1 || x==4)? 3:(x==6)? 1:4))
#define NOTON 5
#define NOTOND 25
#define NOTONT 125
#define EXPLICIT_MAX_LEN 4
#define SINGLE_MIN_LEN 0
#define SINGLE_MAX_LEN 20
#define HAIRPIN_MAX_LEN 30
#define BULGE_MAX_LEN SINGLE_MAX_LEN
#define INTERNAL_MAX_LEN SINGLE_MAX_LEN
#define SYMMETRIC_MAX_LEN 15
#define ASYMMETRY_MAX_LEN 28
#define SPECIAL_HAIRPIN_SCORE_BASELINE -10000

bool _allowed_pairs[NOTON][NOTON];


#define MAXLOOP 30

#define GET_ACGU(x) ((x==1? 'A' : (x==2? 'C' : (x==3? 'G' : (x==4?'U': 'X')))))

#define GET_ACGU_NUC(x) ((x=='A'? 1 : (x=='C'? 2 : (x=='G'? 3 : (x=='U'?4: 0)))))

#define HAIRPINTYPE(x) ((x==5?0 : (x==6?1 : (x==8?2 : 3))))

void initialize()
{
    memset(_allowed_pairs, 0, NOTON * NOTON * sizeof(bool));
    _allowed_pairs[GET_ACGU_NUC('A')][GET_ACGU_NUC('U')] = true;
    _allowed_pairs[GET_ACGU_NUC('U')][GET_ACGU_NUC('A')] = true;
    _allowed_pairs[GET_ACGU_NUC('C')][GET_ACGU_NUC('G')] = true;
    _allowed_pairs[GET_ACGU_NUC('G')][GET_ACGU_NUC('C')] = true;
    _allowed_pairs[GET_ACGU_NUC('G')][GET_ACGU_NUC('U')] = true;
    _allowed_pairs[GET_ACGU_NUC('U')][GET_ACGU_NUC('G')] = true;
}

inline int MIN2(int a, int b) {if (a <= b)return a;else return b;}
inline int MAX2(int a, int b) {if (a >= b)return a;else return b;}

inline int check_special_hairpin(std::string& subseq, int8_t hairpin_type) {
    char *ts;
    const char* tl = subseq.c_str();
    switch(hairpin_type){
        case 0:
            if ((ts=strstr(Triloops, tl))) 
                return -Triloop37[(ts - Triloops)/6];
            break;
        case 1:
            if ((ts=strstr(Tetraloops, tl))) 
                return -Tetraloop37[(ts - Tetraloops)/7];
            break;
        case 2:
            if ((ts=strstr(Hexaloops, tl))) 
                return -Hexaloop37[(ts - Hexaloops)/9];
            break;
        default:
            printf("wrong special hairpin type!\n"); fflush(stdout);
            assert(false);
    }
    return SPECIAL_HAIRPIN_SCORE_BASELINE;
}


inline void v_init_tetra_hex_tri(std::string& seq, int seq_length, std::vector<int>& if_tetraloops, std::vector<int>& if_hexaloops, std::vector<int>& if_triloops) {

    // TetraLoops
    if_tetraloops.resize(seq_length-5<0?0:seq_length-5, -1);
    for (int i = 0; i < seq_length-5; ++i) {
        if (!(seq[i] == 'C' && seq[i+5] == 'G'))
            continue;
        char *ts;
        if ((ts=strstr(Tetraloops, seq.substr(i,6).c_str())))
            if_tetraloops[i] = (ts - Tetraloops)/7;
    }

    // Triloops
    if_triloops.resize(seq_length-4<0?0:seq_length-4, -1);
    for (int i = 0; i < seq_length-4; ++i) {
        if (!((seq[i] == 'C' && seq[i+4] == 'G') || (seq[i] == 'G' && seq[i+4] == 'C')))
            continue;
        char *ts;
        if ((ts=strstr(Triloops, seq.substr(i,5).c_str())))
            if_triloops[i] = (ts - Triloops)/6;
    }

    // Hexaloops
    if_hexaloops.resize(seq_length-7<0?0:seq_length-7, -1);
    for (int i = 0; i < seq_length-7; ++i) {
        if (!(seq[i] == 'A' && seq[i+7] == 'U'))
            continue;
        char *ts;
        if ((ts=strstr(Hexaloops, seq.substr(i,8).c_str())))
            if_hexaloops[i] = (ts - Hexaloops)/9;
    }
    return;
}

inline int v_score_hairpin(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int tetra_hex_tri_index = -1) {
    int size = j-i-1;
    int type = NUC_TO_PAIR(nuci, nucj);

    int energy;

    if(size <= 30)
        energy = hairpin37[size];
    else
        energy = hairpin37[30] + (int)(lxc37*log((size)/30.));

    if(size < 3) return energy;
    #ifdef SPECIAL_HP
    // if(special_hp){
        if (size == 4 && tetra_hex_tri_index > -1)
            return Tetraloop37[tetra_hex_tri_index];
        else if (size == 6 && tetra_hex_tri_index > -1)
            return Hexaloop37[tetra_hex_tri_index];
        else if (size == 3) {
            if (tetra_hex_tri_index > -1)
                return Triloop37[tetra_hex_tri_index];
            return (energy + (type>2 ? TerminalAU37 : 0));
        }
    // }
#endif

    energy += mismatchH37[type][nuci1][nucj_1];

    return energy;
}

inline int v_score_single(int i, int j, int p, int q,
                        int nuci, int nuci1, int nucj_1, int nucj,
                        int nucp_1, int nucp, int nucq, int nucq1){
    int type = NUC_TO_PAIR(nuci, nucj);
    int type_2 = NUC_TO_PAIR(nucq, nucp);
    int n1 = p-i-1;
    int n2 = j-q-1;
    int nl, ns, u, energy;
    energy = 0;

    if (n1>n2) { nl=n1; ns=n2;}
    else {nl=n2; ns=n1;}

    if (nl == 0)
        return stack37[type][type_2];  /* stack */

    if (ns==0) {                      /* bulge */
        energy = (nl<=MAXLOOP)?bulge37[nl]:
      (bulge37[30]+(int)(lxc37*log(nl/30.)));
    if (nl==1) energy += stack37[type][type_2];
    else {
      if (type>2) energy += TerminalAU37;
      if (type_2>2) energy += TerminalAU37;
    }
    return energy;
  }
  else {                            /* interior loop */
    if (ns==1) {
      if (nl==1)                    /* 1x1 loop */
        return int11_37[type][type_2][nuci1][nucj_1];
      if (nl==2) {                  /* 2x1 loop */
        if (n1==1)
          energy = int21_37[type][type_2][nuci1][nucq1][nucj_1];
        else
          energy = int21_37[type_2][type][nucq1][nuci1][nucp_1];
        return energy;
      }
      else {  /* 1xn loop */
        energy = (nl+1<=MAXLOOP)?(internal_loop37[nl+1]) : (internal_loop37[30]+(int)(lxc37*log((nl+1)/30.)));
        energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);
        energy += mismatch1nI37[type][nuci1][nucj_1] + mismatch1nI37[type_2][nucq1][nucp_1];
        return energy;
      }
    }
    else if (ns==2) {
      if(nl==2)      {              /* 2x2 loop */
        return int22_37[type][type_2][nuci1][nucp_1][nucq1][nucj_1];}
      else if (nl==3){              /* 2x3 loop */
        energy = internal_loop37[5]+ninio37;
        energy += mismatch23I37[type][nuci1][nucj_1] + mismatch23I37[type_2][nucq1][nucp_1];
        return energy;
      }

    }
    { /* generic interior loop (no else here!)*/
      u = nl + ns;
      energy = (u <= MAXLOOP) ? (internal_loop37[u]) : (internal_loop37[30]+(int)(lxc37*log((u)/30.)));

      energy += MIN2(MAX_NINIO, (nl-ns)*ninio37);

      energy += mismatchI37[type][nuci1][nucj_1] + mismatchI37[type_2][nucq1][nucp_1];
    }
  }
  return energy;
}

// multi_loop
inline int E_MLstem(int type, int si1, int sj1) {
    int energy = 0;

    if(si1 >= 0 && sj1 >= 0){
        energy += mismatchM37[type][si1][sj1];
    }
    else if (si1 >= 0){
        energy += dangle5_37[type][si1];
    }
    else if (sj1 >= 0){
        energy += dangle3_37[type][sj1];
    }

    if(type > 2) {
        energy += TerminalAU37;
    }

    energy += ML_intern37;

    return energy;
}

// hzhang: added for mRNA design
inline int E_MLstem_without_Dangle(int type, int si1, int sj1) {
    int energy = 0;

    if(type > 2) {
        energy += TerminalAU37;
    }

    energy += ML_intern37;

    return energy;
}

inline int v_score_M1(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
    int p = i;
    int q = k;
    int tt = NUC_TO_PAIR(nuci, nuck);

    return E_MLstem(tt, nuci_1, nuck1);

}

inline int v_score_M1_without_dangle(int i, int j, int k, int nuci_1, int nuci, int nuck, int nuck1, int len) {
    int tt = NUC_TO_PAIR(nuci, nuck);

    return E_MLstem_without_Dangle(tt, -1, -1);

}

inline int v_score_multi_unpaired(int i, int j) {
    return 0;
}

inline int v_score_multi(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    int tt = NUC_TO_PAIR(nucj, nuci);

    return E_MLstem(tt, nucj_1, nuci1) + ML_closing37;
}

inline int v_score_multi_without_dangle(int i, int j, int nuci, int nuci1, int nucj_1, int nucj, int len) {
    int tt = NUC_TO_PAIR(nucj, nuci);

    return E_MLstem_without_Dangle(tt, -1, -1) + ML_closing37;
}

// exterior_loop
inline int v_score_external_paired(int i, int j, int nuci_1, int nuci, int nucj, int nucj1, int len) {
    int type = NUC_TO_PAIR(nuci, nucj);
    int energy = 0;

    if(nuci_1 >= 0 && nucj1 >= 0){
        energy += mismatchExt37[type][nuci_1][nucj1];
    }
    else if (nuci_1 >= 0){
        energy += dangle5_37[type][nuci_1];
    }
    else if (nucj1 >= 0){
        energy += dangle3_37[type][nucj1];
    }

    if(type > 2)
        energy += TerminalAU37;
  return energy;
}

inline int v_score_external_paired_without_dangle(int i, int j, int nuci, int nucj, int len) {
    int type = NUC_TO_PAIR(nuci, nucj);
    int energy = 0;

    if(type > 2)
        energy += TerminalAU37;
  return energy;
}

inline int v_score_external_unpaired(int i, int j) {
    return 0;
}

#endif //FASTCKY_UTILITY_V_H

