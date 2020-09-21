// Copyright (C) 2015-2020 Bhalchandra Digambar Thatte

// This file is a part of the project ASSIGN - a program to allocate disciplinas
// (subjects) to professors. ASSIGN project is released under the GPL V3 license.

// Some problems to be solved

// 1. There are turmas (sections) with different codes with schedule conflict,
// e.g., N and N1. They cannot be assigned the same code because some people may
// have preference for N over N1. Chanding N1 to Y can potentially lead to
// solutions that are not feasible; e.g., the program may assign 38N and 38Y to
// a person, since it will not detect schedule conflict between 38N and 38Y. Of
// course, this will be rare since 38N and 38Y will never be a person's
// preference. One possible solution is to have groups of conflicting codes
// instead of assuming that codes are mutually non-conflicting.

// 2. function cost_map_add: The cost of i-th preference for 0 <= i <= 5 is x +
// w_i y, where x and y depend as described next on the min hours the prof is
// teaching. The idea is to value more the preferences of people teaching
// more. Perhaps we should simply use weights since people teaching less may
// have already been allocated other turmas (sections).

#include <algorithm>
#include <array>
#include <assert.h>
#include <cstdlib> // std::rand, std::srand
#include <ctime> // std::time
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock

#define COST_ID 8 // cost of teaching different subjects
#define COST_LOW 100 // cost of not teaching enough hours
#define COST_HIGH 100 // cost of teaching more hours
#define COST_PREF_HOURLY 60
#define COST_PREF_MIN 200 // min cost of violating preferences
#define COST_MAX 1000 // schedule conflict

// Currently preferences are numbered 0,...,5. If a given combination is
// not a preference of a prof, assign the allocation a preference no. 6
#define MAX_PREF 6

// --------------------------------------------------------------------
// Brown 2013, Random number generation in C++11.
std::default_random_engine & global_urng(){
  static std::default_random_engine u{};
  return u;
}
void randomize(int s = 0){
  static std::random_device rd;
  if (s) global_urng().seed(s);
  else global_urng().seed(rd());
}
// random double from the interval [lo,hi)
double random_real(double lo, double hi){
  static std::uniform_real_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(global_urng(),parm_t{lo,hi});
}
// random int from the interval [lo, hi)
int random_int(int lo, int hi){
  static std::uniform_int_distribution<> d{};
  using parm_t = decltype(d)::param_type;
  return d(global_urng(),parm_t{lo,hi-1});
}
// --------------------------------------------------------------------

struct Assignment;
struct Prof;
struct Disc;

// Global variables
uint8_t ND = 0; // number of disciplinas
uint8_t NP = 0; // number of profs
std::vector<Disc> DV; // disciplinas
std::map<std::string,uint8_t> DM; // <disciplina, index in D>
std::vector<Prof> PV;
std::map<Assignment,unsigned int> cost_map;

#ifdef Assignment_ARR
// Assignment is implemented as an array of 4 ints
struct Assignment {
  std::array<uint8_t,4> arr;
public:
  Assignment(){arr = {};}
  Assignment(std::array<uint8_t,4>& v){
    arr.swap(v);
    std::sort(arr.begin(),arr.begin()+3);
  }
  void operator=(Assignment b){
    arr.swap(b.arr);
    std::sort(arr.begin(),arr.begin()+3);
  }
  uint8_t operator[](unsigned int i){return arr[i];}
  bool operator!=(const Assignment& v){return (arr != v.arr);}
  std::string print();
  unsigned int hours(); // total number of hours of this assignment
  uint32_t ass2uint32() const{
    return (((unsigned int)arr[0] << 24) |
	    ((unsigned int)arr[1] << 16) |
	    ((unsigned int)arr[2] << 8) |
	    ((unsigned int)arr[3]));
  }
};

inline bool operator==(const Assignment& lhs, const Assignment& rhs){
  return (lhs.arr == rhs.arr);
}
inline bool operator<(const Assignment& lhs, const Assignment& rhs){
  return (lhs.ass2uint32() < rhs.ass2uint32());
}

#else // Assignment implemented as a single uint32_t

struct Assignment {
  uint32_t v;
public:
  Assignment(){};
  Assignment(std::array<uint8_t,4>& u){
    std::sort(u.begin(),u.begin()+3);
    v = (((uint32_t) u[0]) << 24) |
      (((uint32_t) u[1]) << 16) |
      (((uint32_t) u[2]) << 8) |
      ((uint32_t) u[3]);
  };
  void operator=(Assignment b){
    v = b.v;
  }
  uint8_t operator[](unsigned int i){
    return ((v >> (8*(3-i))) & 255);
  }
  bool operator!=(const Assignment& u){return (v != u.v);}
  std::string print();
  unsigned int hours(); // total number of hours of this assignment
};

inline bool operator==(const Assignment& lhs, const Assignment& rhs){
  return (lhs.v == rhs.v);
}
inline bool operator<(const Assignment& lhs, const Assignment& rhs){
  return (lhs.v < rhs.v);
}
#endif

// Professor
struct Prof {
  std::string name; //name
  unsigned int hmin; // min. number of hours to teach
  unsigned int hmax; // max. number of hours to teach
  std::vector<Assignment> prefs; // list of preferences
  Assignment current_ass; // current assignment
  unsigned int current_cost=0; // cost of assignment; 0 means nothing assigned yet
  unsigned int min_pref_cost = COST_PREF_MIN;
  bool active0 = true;
  bool active1 = true;
  uint8_t snp; // short name
public:
  void set_ass(Assignment ass, unsigned int c = 0){
    current_ass = ass;
    current_cost = c;
    if (c == 0) {
      auto it = cost_map.find(ass);
      if(it != cost_map.end()){
        current_cost = it->second;
      }
      else
	/// should not be here
	current_cost = COST_MAX;
    }
  };
  std::string print();
  int ass_pref();
  int ass_pref(Assignment& ass);
};

// Disciplina (subject)
struct Disc {
  unsigned int id;
  char code; // A,B,C,D,E,F,M,N,Z: denotes timetable slot
  unsigned int nt; // number of turmas (sections)
  unsigned int nh; // hours per week
  unsigned int free; // not assigned yet
  uint8_t snd; // short name of a disciplina
};


#ifdef Assignment_ARR
// Assignment is implemented as an array of 4 ints
std::string Assignment::print(){
  std::ostringstream ost;
  for(size_t i = 0; i < 3; i++){
    uint8_t j = arr[i];
    if (j){
	  if(i < 2)
		ost << DV[j].id << DV[j].code << " ";
	  else if(2==i)
		ost << DV[j].id << DV[j].code;
    }
  }
  return ost.str();
}

unsigned int Assignment::hours(){
  unsigned int h = 0;
  for(size_t i = 0; i < 3; ++i)
    h+= DV[arr[i]].nh;
  return h;
}

#else
// Assignment is implemented as a single uint32_t
std::string Assignment::print(){
  std::ostringstream ost;
  for(size_t i = 3; i > 0; i--){
    unsigned int j = (v >> 8*i) & 255;
    if (j){
	  if(i>1)
		ost << DV[j].id << DV[j].code << " ";
	  else if(1==i)
		ost << DV[j].id << DV[j].code;
    }
  }
  return ost.str();
}
unsigned int Assignment::hours(){
  unsigned int h = 0;
  for(size_t i = 3; i > 0; i--){
    unsigned int j = (v >> 8*i) & 255;
    h+= DV[j].nh;
  }
  return h;
}
#endif

std::string Prof::print(){
  std::ostringstream ost;
  ost << name << ", " << hmin << ", "  << hmax << ", ";
  for(auto it = prefs.begin(); it != prefs.end(); ++it){
	ost << (*it).print() << ", ";
  }
  return ost.str();
}

int Prof::ass_pref(){
  int i = MAX_PREF;
  for(auto it = prefs.begin(); it != prefs.end(); ++it){
	if (current_ass == (*it)){
	  i = it-prefs.begin();
	  break;
	}
  }
  return i;
}

int Prof::ass_pref(Assignment& ass){
  int i = MAX_PREF;
  for(auto it = prefs.begin(); it != prefs.end(); ++it){
	if (ass == (*it)){
	  i = it-prefs.begin();
	  break;
	}
  }
  return i;
}

// Total cost of an assignment for all professors
unsigned int total_cost(){
  unsigned int total_cost = 0;
  for(auto it = PV.begin(); it != PV.end(); ++it){
    total_cost += (*it).current_cost;
  }
  return total_cost;
}

// Average assignment preference
double average_pref(){
  unsigned int total_pref = 0;
  unsigned int unmatched = 0;
  for(auto it = PV.begin(); it != PV.end(); ++it){
    int i = (*it).ass_pref();
    if(i == MAX_PREF)
      unmatched++;
    else
      total_pref += i;
  }
  return (double)total_pref/((double)(NP-unmatched));
}

void print_solution(std::string msg){
  // output allocation
  auto tc = total_cost();
  std::cout << "\n";
  std::cout << msg << "\n";
  std::cout
    << std::setw(32) << "Name" << ","
    << std::setw(8) << "min_h" << ","
    << std::setw(8) << "max_h" << ","
    << std::setw(16) << "Allocation" << ","
    << std::setw(6) << "Pref" << ","
    << std::setw(6) << "Cost"
    << std::endl;

  for(auto it = PV.begin(); it != PV.end(); ++it){
    std::cout
      << std::setw(32) << (*it).name << ","
      << std::setw(8) << (*it).hmin << ","
      << std::setw(8) << (*it).hmax << ","
      << std::setw(16) << (*it).current_ass.print() << ","
      << std::setw(6) << (*it).ass_pref() << ","
      << std::setw(6) << (*it).current_cost
      << std::endl;
  }
  std::cout << "Total cost = " << tc << "\n";
  std::cout << "Average preference = " << average_pref() << "\n";
  return;
}

bool schedule_conflict(Assignment& a){
  uint8_t i = a[0];
  uint8_t j = a[1];
  uint8_t k = a[2];

  if (i == 0 && j == 0)	return false;
  else if ((i == 0) && (j != 0) && (k != 0)){
    if (DV[j].code == DV[k].code) return true;
  }	else if ((DV[i].code == DV[j].code)	|| (DV[i].code == DV[k].code) || (DV[j].code == DV[k].code))
    return true;
  return false;
}

void cost_map_add(Assignment& c){
  // Return if already in the map.
  auto it = cost_map.find(c);
  if(it != cost_map.end()){
    return;
  }

  if(schedule_conflict(c)){
    cost_map[c] = COST_MAX;
    return;
  }

  unsigned int total_cost = 0;
  uint8_t pi = c[3];
  auto hmin = PV[pi].hmin;

  // If Assignment c is preferred by prof pi, set preference cost and
  // return; cost of i-th preference for 0 <= i <= 5 is x + w_i y, where
  // x and y depend as described next on the min hours the prof is
  // teaching. The idea is to value more the preferences of people
  // teaching more. (See TODO item 2.)

  int pref = PV[pi].ass_pref(c);
  if(pref < MAX_PREF) { // pref is one of the preferences
    // cost of i-th preference is (x+iy)
    unsigned int x = 0;
    unsigned int y = 1;
    if(hmin >= 8) {x = 3; y = 2;}
    if(hmin >= 12){x = 6; y = 3;}
    cost_map[c] = x + pref*y;
    return;
  } // else (pref == MAX_PREF): add cost of non-preference

  auto h = c.hours();
  unsigned int pref_cost;
  if(PV[pi].min_pref_cost == 0) pref_cost = 0;
  else {
    pref_cost = h * COST_PREF_HOURLY;
    pref_cost = (pref_cost < PV[pi].min_pref_cost) ? PV[pi].min_pref_cost : pref_cost;
  }
  total_cost += pref_cost;

  // Check if h within [hmin,hmax] range

  if(h < PV[pi].hmin){
    total_cost += COST_LOW;
  } else if (h > PV[pi].hmax) {
    total_cost += COST_HIGH * pow(2,h-PV[pi].hmax);
  }

  uint8_t i = c[0];
  uint8_t j = c[1];
  uint8_t k = c[2];

  if (i == 0 && j == 0){
    // assignment has at most one non-zero turma
    cost_map[c] = total_cost;
    return;
  }

  // => assignment has 2 or 3 turmas (sections)
  // Cost of teaching different subjects
  if ((i == 0) && (j != 0) && (k != 0)){ // => assignment has 2 turmas (sections)
    if ((DV[j].id != DV[k].id)){
      // => assignment has different subjects
      total_cost += COST_ID;
    }
  } else { // => assignment of 3 turmas (sections)
    if (((DV[i].id != DV[j].id) || (DV[i].id != DV[k].id) || (DV[j].id != DV[k].id))){
      total_cost += COST_ID;
    }
  }
  cost_map[c] = total_cost;
  return;
}

// Compute costs of all possible individual assignments, i.e., roughly
// NP * binom{ND}{3} and store them in a map<Assignment,cost> for quick look
// up during optimisation.
void preprocess(){
  for(auto pit = PV.begin(); pit != PV.end(); ++pit){
    for(size_t i = 0; i < ND; ++i){
      for(size_t j = i; j < ND; ++j){
        for(size_t k = j; k < ND; ++k){
          std::array<uint8_t,4> z = {DV[i].snd,DV[j].snd,DV[k].snd,(*pit).snp};
          Assignment zc(z);
          cost_map_add(zc);
        } // k
      } // j
    } // i
  } // pit
}

// Initial solution. Randomises PV, assigns turmas (sections) to profs in a
// greedy manner, and sets current_ass and current_cost for each prof.
void init() {
  // std::random_shuffle(PV.begin(), PV.end());
  std::shuffle(PV.begin(), PV.end(), global_urng());
  // reset
  for(auto pit = PV.begin(); pit != PV.end(); ++pit){
    std::array<uint8_t,4> x = {0,0,0,(*pit).snp};
    Assignment xa(x);
    (*pit).set_ass(xa,cost_map[xa]);
    (*pit).active0 = true;
    (*pit).active1 = true;
  }
  for(auto dit = DV.begin(); dit != DV.end(); ++dit){
    (*dit).free = (*dit).nt;
  }
  std::vector<Prof>::iterator pit = PV.begin();
  std::vector<Disc>::iterator dit = DV.begin();
  ++dit; // skip the null discipline

  while (pit != PV.end() && dit != DV.end()) {
    unsigned int h = 0;
    std::array<uint8_t,4> x = {};
    size_t i = 0;
    while (h+4 <= (*pit).hmax && i < 3) {
      if ((*dit).free > 0) {
	x[i++] = dit-DV.begin();
	h+= (*dit).nh;
	--((*dit).free);
      }
      if((*dit).free == 0){
	++dit;
	if(dit == DV.end()) break;
      }
    } // => (h+4 > (*pit).max) or i = 3 or dit == DV.end()

    x[3] = (*pit).snp;
    Assignment xa(x);
    (*pit).set_ass(xa,cost_map[xa]);
    ++pit;
  }
  if(pit != PV.end() && dit == DV.end()){
    while (pit != PV.end()){
      std::array<uint8_t,4> x = {0,0,0,(*pit).snp};
      Assignment xa(x);
      (*pit).set_ass(xa,cost_map[xa]);
      ++pit;
    }
  }
  if(pit == PV.end() && dit != DV.end()){
    std::cout << "not enough profs" << std::endl;
  }
}

// Try swapping all possible combinations of turmas (sections) in the current
// assignments of 4 profs. Returns delta (change in total cost) (0 if swap not
// made, -ve if swap is made). accept_prob is the probability of accepting a bad
// solution (+ve in case of simulated annealing, 0 otherwise).
int swap4(unsigned int pi, unsigned int pj, unsigned int pk, unsigned int pl,
          double accept_prob = 0){
  Assignment pi_best_ass = PV[pi].current_ass; // current assignment to pi
  Assignment pj_best_ass = PV[pj].current_ass; // current assignment to pj
  Assignment pk_best_ass = PV[pk].current_ass; // current assignment to pk
  Assignment pl_best_ass = PV[pl].current_ass; // current assignment to pl
  auto pi_min_cost = PV[pi].current_cost; // current cost to pi
  auto pj_min_cost = PV[pj].current_cost; // current cost to pj
  auto pk_min_cost = PV[pk].current_cost; // current cost to pk
  auto pl_min_cost = PV[pl].current_cost; // current cost to pl

  std::array<uint8_t,12> t =
    {
      PV[pi].current_ass[0],PV[pi].current_ass[1],PV[pi].current_ass[2],
      PV[pj].current_ass[0],PV[pj].current_ass[1],PV[pj].current_ass[2],
      PV[pk].current_ass[0],PV[pk].current_ass[1],PV[pk].current_ass[2],
      PV[pl].current_ass[0],PV[pl].current_ass[1],PV[pl].current_ass[2]
    };
  std::sort(t.begin(),t.end());
  uint8_t m[12] = {0}; // mask

  bool commit = false;
  int delta_min = 0;
  for(size_t i = 0; i < 10; ++i){
    for(size_t j = i+1; j < 11; ++j){
      for(size_t k = j+1; k < 12; ++k){

        std::array<uint8_t,4> x = {t[i],t[j],t[k],PV[pi].snp};
        Assignment pi_ass(x);
        if(schedule_conflict(pi_ass)) continue;
        auto pi_cost = cost_map[pi_ass];
        if(pi_cost >= pi_min_cost + COST_PREF_MIN) continue;
        m[i] = m[j] = m[k] = 1;

        for(size_t i1 = 0; i1 < 10; ++i1){
          if(m[i1]) continue;
          for(size_t j1 = i1+1; j1 < 11; ++j1){
            if(m[j1]) continue;
            for(size_t k1 = j1+1; k1 < 12; ++k1){
              if(m[k1]) continue;

              std::array<uint8_t,4> y = {t[i1],t[j1],t[k1],PV[pj].snp};
              Assignment pj_ass(y);
              if(schedule_conflict(pj_ass)) continue;
              auto pj_cost = cost_map[pj_ass];
              if(pj_cost >= pj_min_cost + COST_PREF_MIN) continue;
              m[i1] = m[j1] = m[k1] = 2;

              for(size_t i2 = 0; i2 < 10; ++i2){
                if(m[i2]) continue;
                for(size_t j2 = i2+1; j2 < 11; ++j2){
                  if(m[j2]) continue;
                  for(size_t k2 = j2+1; k2 < 12; ++k2){
                    if(m[k2]) continue;
                    std::array<uint8_t,4> z = {t[i2],t[j2],t[k2],PV[pk].snp};
                    Assignment pk_ass(z);
                    if(schedule_conflict(pk_ass)) continue;
                    auto pk_cost = cost_map[pk_ass];
                    if(pk_cost >= pk_min_cost + COST_PREF_MIN) continue;
                    m[i2] = m[j2] = m[k2] = 3;

                    std::array<uint8_t,4> w = {0,0,0,PV[pl].snp};
                    size_t wi = 0;
                    for(size_t l = 0; l < 12; l++){
                      if(m[l] == 0) w[wi++] = t[l];
                    }
                    Assignment pl_ass(w);
                    if(schedule_conflict(pl_ass)) {
                      m[k2]=0;
                      continue;
                    }
                    auto pl_cost = cost_map[pl_ass];

                    int delta = pi_cost - PV[pi].current_cost + pj_cost - PV[pj].current_cost
                      + pk_cost - PV[pk].current_cost + pl_cost - PV[pl].current_cost;

                    // accept a good solution always, and a bad solution with
                    // probability accept_prob
                    auto r = 0.0;
                    if(accept_prob > 0.0)
                      r = random_real(0.0,1.0);
                    if ((delta < delta_min) || (r < accept_prob)){
                      pi_best_ass = pi_ass;
                      pj_best_ass = pj_ass;
                      pk_best_ass = pk_ass;
                      pl_best_ass = pl_ass;
                      pi_min_cost = pi_cost;
                      pj_min_cost = pj_cost;
                      pk_min_cost = pk_cost;
                      pl_min_cost = pl_cost;
                      delta_min = delta;
                      commit = true;
                    }
                    m[k2]=0;
                  } // k2
                  m[j2]=0;
                } // j2
                m[i2]=0;
              } // i2
              m[k1]=0;
            } // k1
            m[j1]=0;
          } // j1
          m[i1]=0;
        } // i1
        m[k]=0;
      } //k
      m[j]=0;
    }//j
    m[i]=0;
  }//i

  if (commit) {
    if(pi_best_ass != PV[pi].current_ass){
      PV[pi].set_ass(pi_best_ass, pi_min_cost);
      PV[pi].active1 = true;
    }
    if(pj_best_ass != PV[pj].current_ass){
      PV[pj].set_ass(pj_best_ass, pj_min_cost);
      PV[pj].active1 = true;
    }
    if(pk_best_ass != PV[pk].current_ass){
      PV[pk].set_ass(pk_best_ass, pk_min_cost);
      PV[pk].active1 = true;
    }
    if(pl_best_ass != PV[pl].current_ass){
      PV[pl].set_ass(pl_best_ass, pl_min_cost);
      PV[pl].active1 = true;
    }
  }
  if (commit) return -1;
  else return 1;
}

// Try swapping all possible combinations of turmas (sections) in the current
// assignments of 3 profs. Returns delta (change in total cost) (0 if swap not
// made, -ve if swap is made). accept_prob is the probability of accepting a bad
// solution (+ve in case of simulated annealing, 0 otherwise).
int swap3(unsigned int pi, unsigned int pj, unsigned int pk,
          double accept_prob=0){
  Assignment pi_best_ass = PV[pi].current_ass; // current assignment to pi
  Assignment pj_best_ass = PV[pj].current_ass; // current assignment to pj
  Assignment pk_best_ass = PV[pk].current_ass; // current assignment to pk
  auto pi_min_cost = PV[pi].current_cost; // current cost to pi
  auto pj_min_cost = PV[pj].current_cost; // current cost to pj
  auto pk_min_cost = PV[pk].current_cost; // current cost to pk

  std::array<uint8_t,9> t = {PV[pi].current_ass[0],PV[pi].current_ass[1],PV[pi].current_ass[2],
                             PV[pj].current_ass[0],PV[pj].current_ass[1],PV[pj].current_ass[2],
                             PV[pk].current_ass[0],PV[pk].current_ass[1],PV[pk].current_ass[2]};
  std::sort(t.begin(),t.end());
  uint8_t m[9] = {0}; // mask

  bool commit = false;
  int delta_min = 0;
  for(size_t i = 0; i < 7; ++i){
    for(size_t j = i+1; j < 8; ++j){
      for(size_t k = j+1; k < 9; ++k){

        std::array<uint8_t,4> x = {t[i],t[j],t[k],PV[pi].snp};
        Assignment pi_ass(x);
        if(schedule_conflict(pi_ass)) continue;
        auto pi_cost = cost_map[pi_ass];
        if(pi_cost > pi_min_cost + COST_PREF_MIN) continue;
        m[i] = m[j] = m[k] = 1;

        for(size_t i1 = 0; i1 < 7; ++i1){
          if(m[i1]) continue;
          for(size_t j1 = i1+1; j1 < 8; ++j1){
            if(m[j1]) continue;
            for(size_t k1 = j1+1; k1 < 9; ++k1){
              if(m[k1]) continue;

              std::array<uint8_t,4> y = {t[i1],t[j1],t[k1],PV[pj].snp};
              Assignment pj_ass(y);
              if(schedule_conflict(pj_ass)) continue;
              auto pj_cost = cost_map[pj_ass];
              if(pj_cost > pj_min_cost + COST_PREF_MIN)  continue;
              m[i1] = m[j1] = m[k1] = 2;

              std::array<uint8_t,4> z = {0,0,0,PV[pk].snp};
              size_t zi = 0;
              for(size_t l = 0; l < 9; l++){
                if(m[l] == 0) z[zi++] = t[l];
              }
              Assignment pk_ass(z);
              if(schedule_conflict(pk_ass)){
                m[k1] = 0;
                continue;
              }
              auto pk_cost = cost_map[pk_ass];

              int delta = pi_cost - PV[pi].current_cost + pj_cost - PV[pj].current_cost
                + pk_cost - PV[pk].current_cost;

              // accept a good solution always, and a bad solution with
              // probability accept_prob
              auto r = 0.0;
              if(accept_prob > 0.0)
                r = random_real(0.0,1.0);
              if ((delta < delta_min) || (r < accept_prob)){
                pi_min_cost = pi_cost;
                pi_best_ass = pi_ass;
                pj_best_ass = pj_ass;
                pj_min_cost = pj_cost;
                pk_best_ass = pk_ass;
                pk_min_cost = pk_cost;
                delta_min = delta;
                commit = true;
              }
              m[k1]=0;
            } // k1
            m[j1]=0;
          } // j1
          m[i1]=0;
        } // i1
        m[k]=0;
      } // k
      m[j]=0;
    } // j
    m[i]=0;
  } // i

  if (commit) {
    if(pi_best_ass != PV[pi].current_ass){
      PV[pi].set_ass(pi_best_ass, pi_min_cost);
      PV[pi].active1 = true;
    }
    if(pj_best_ass != PV[pj].current_ass){
      PV[pj].set_ass(pj_best_ass, pj_min_cost);
      PV[pj].active1 = true;
    }
    if(pk_best_ass != PV[pk].current_ass){
      PV[pk].set_ass(pk_best_ass, pk_min_cost);
      PV[pk].active1 = true;
    }
  }
  if (commit) return -1;
  else return 1;
}

// Try swapping all possible combinations of turmas (sections) in the current
// assignments of 2 profs. Returns delta (change in total cost) (0 if swap not
// made, -ve if swap is made). accept_prob is the probability of accepting a bad
// solution (+ve in case of simulated annealing, 0 otherwise).
int swap2(unsigned int pi, unsigned int pj, double accept_prob=0){
  bool commit = false;
  int delta_min = 0;
  Assignment pi_best_ass = PV[pi].current_ass; // current assignment to pi
  Assignment pj_best_ass = PV[pj].current_ass; // current assignment to pj
  auto pi_min_cost = PV[pi].current_cost; // current cost to pi
  auto pj_min_cost = PV[pj].current_cost; // current cost to pj

  std::array<uint8_t,6> t = {PV[pi].current_ass[0],PV[pi].current_ass[1],PV[pi].current_ass[2],
			     PV[pj].current_ass[0],PV[pj].current_ass[1],PV[pj].current_ass[2]};
  std::sort(t.begin(),t.end());
  // try all combinations of possible swaps
  for(size_t i = 0; i < 4; ++i){
    for(size_t j = i+1; j < 5; ++j){
      for(size_t k = j+1; k < 6; ++k){

        std::array<uint8_t,4> x = {t[i],t[j],t[k],PV[pi].snp};
        Assignment pi_ass(x);
        if(schedule_conflict(pi_ass)) continue;
        auto pi_cost = cost_map[pi_ass];
        if(pi_cost > pi_min_cost + COST_PREF_MIN) continue;

        std::array<uint8_t,4> y = {0,0,0,PV[pj].snp};
        size_t yi=0;
        for(size_t m = 0; m < 6; m++){
          if ((m != i) && (m != j) && (m != k)) {
            y[yi++] = t[m];
          }
        }
        Assignment pj_ass(y);
        if(schedule_conflict(pj_ass)) continue;
        auto pj_cost = cost_map[pj_ass];
        if(pj_cost > pj_min_cost + COST_PREF_MIN)  continue;

        int delta = pi_cost - PV[pi].current_cost + pj_cost - PV[pj].current_cost;

        // accept a good solution always, and a bad solution with
        // probability accept_prob
        auto r = 0.0;
        if(accept_prob > 0.0)
          r = random_real(0.0,1.0);
        if ((delta < delta_min) || (r < accept_prob)){
          pi_min_cost = pi_cost;
          pi_best_ass = pi_ass;
          pj_best_ass = pj_ass;
          pj_min_cost = pj_cost;
          delta_min = delta;
          commit = true;
        }
      } // k
    } // j
  } // i

  if (commit) {
    if(pi_best_ass != PV[pi].current_ass){
      PV[pi].set_ass(pi_best_ass, pi_min_cost);
      PV[pi].active1 = true;
    }
    if(pj_best_ass != PV[pj].current_ass){
      PV[pj].set_ass(pj_best_ass, pj_min_cost);
      PV[pj].active1 = true;
    }
  }
  if (commit) return -1;
  else return 1;
}

// Hill climbing with swap2 - calls swap2 for all pairs of profs
// repeatedly until a local minimum is reached.
void opt2(double accept_prob=0){
  // reset active0, active1 flags
  for(auto it = PV.begin(); it != PV.end(); ++it){
    (*it).active0 = true;
    (*it).active1 = true;
  }
  bool stable = false;
  while (!stable) {
    stable = true;
    // For each pair of profs, check if swapping some disciplinas
    // between them can lower cost
    for(size_t i = 0; i < NP; ++i){
      for(size_t j = i+1; j < NP; ++j){
	// Call swap fori, j if at least one of them is active (changed
	// in the previous iteration); otherwise no need to compare i,j
	// again.
	if(PV[i].active0 || PV[j].active0){
	  // Set stable = false if swap happens fori,j
	  if(swap2(i,j, accept_prob) < 0) // swap achieved reduction in cost
	    stable = false;
	}
      }
    }
    if (stable) break; // no swap happened in the loops
    else {
      for(auto it = PV.begin(); it != PV.end(); ++it){
	(*it).active0 = (*it).active1;
	(*it).active1 = false;
      }
    }
  }
  assert(stable);
  // print_solution("Algo: opt2");
}

// Hill climbing with swap3 - calls swap3 for all triples of profs
// repeatedly until a local minimum is reached.
void opt3(double accept_prob=0){
  // reset active0, active1 flags
  for(auto it = PV.begin(); it != PV.end(); ++it){
    (*it).active0 = true;
    (*it).active1 = true;
  }
  bool stable = false;
  while (!stable) {
    stable = true;
    // For each triple of profs, check if swapping some disciplinas
    // between them can lower cost
    for(size_t i = 0; i < NP; ++i){
      for(size_t j = i+1; j < NP; ++j){
        for(size_t k = j+1; k < NP; ++k){
          // Call swap fori, j, k if at least one of them is active
          // (changed in the previous iteration); otherwise no need to
          // compare i,j,k again.
          if(PV[i].active0 || PV[j].active0 || PV[k].active0){
            // Set stable = false if swap happens fori,j,k
            if(swap3(i,j,k, accept_prob) < 0) // swap achieved reduction in cost
              stable = false;
          }
        }// k
      }// j
    }// i
    if (stable) break; // no swap happened in the loops
    else {
      for(auto it = PV.begin(); it != PV.end(); ++it){
        (*it).active0 = (*it).active1;
        (*it).active1 = false;
      }
    }

  }
  assert(stable);
  // print_solution("Algo: opt3");
}

// Hill climbing with swap4 - calls swap4 for all quartets of profs
// repeatedly until a local minimum is reached.
void opt4(double accept_prob=0){
  // reset active0, active1 flags
  for(auto it = PV.begin(); it != PV.end(); ++it){
    (*it).active0 = true;
    (*it).active1 = true;
  }
  bool stable = false;
  while (!stable) {
    stable = true;
    // For each quadruple of profs, check if swapping some disciplinas
    // between them can lower cost
    for(size_t i = 0; i < NP; ++i){
      for(size_t j = i+1; j < NP; ++j){
        for(size_t k = j+1; k < NP; ++k){
          for(size_t l = k+1; l < NP; l++){
            // Call swap for i, j, k, l if at least one of them is active
            // (changed in the previous iteration); otherwise no need to
            // compare i,j,k,l again.
            if(PV[i].active0 || PV[j].active0 || PV[k].active0 || PV[l].active0){
              // Set stable = false if swap happens fori,j,k,l
              if(swap4(i,j,k,l, accept_prob) < 0) // swap achieved reduction in cost
                stable = false;
            }
          }//l
        }//k
      }//j
    }//i
    if (stable) break; // no swap happened in the loops
    else {
      for(auto it = PV.begin(); it != PV.end(); ++it){
        (*it).active0 = (*it).active1;
        (*it).active1 = false;
      }
    }
  }
  assert(stable);
  //  print_solution("Algo: opt4");
  return;
}

// Hill climbing with swap3 (mostly) with in between steps of swap4
// until the solution is perturbed. This does not seem to work very
// well.
void opt34(){
  for(auto it = PV.begin(); it != PV.end(); ++it){
    (*it).active0 = true;
    (*it).active1 = true;
  }
  bool stable3 = false;
  bool stable4 = false;

  while ((!stable3) || (!stable4)) {
    stable3 = true;
    // For each quadruple pair of profs, check if swapping some
    // disciplinas between them can lower cost
    for(size_t i = 0; i < NP; ++i){
      for(size_t j = i+1; j < NP; ++j){
        for(size_t k = j+1; k < NP; ++k){
          // Call swap for i, j, k if at least one of them is active
          // (changed in the previous iteration); otherwise no need to
          // compare i,j,k again.
          if(PV[i].active0 || PV[j].active0 || PV[k].active0){
            // Set stable = false if swap happens for i,j,k
            if(swap3(i,j,k) < 0) // swap achieved reduction in cost
              stable3 = false;
          }
        }//k
      }//j
    }//i
    if (!stable3) { // some swap3 steps succeeded
      for(auto it = PV.begin(); it != PV.end(); ++it){
        (*it).active0 = (*it).active1;
        (*it).active1 = false;
      }
    } else {
      print_solution("Algo 3");
      // reset active nodes
      for(auto it = PV.begin(); it != PV.end(); ++it){
        (*it).active0 = true;
        (*it).active1 = true;
      }

      // do swap4 until it becomes unstable again
      bool stable4 = true;
      while (stable4) {
	// foreach quadruple pair of profs, check if swapping some
	// disciplinas between them can lower cost
        for(size_t i = 0; i < NP; ++i){
          for(size_t j = i+1; j < NP; ++j){
            for(size_t k = j+1; k < NP; ++k){
              for(size_t l = k+1; l < NP; l++){
                // Call swap fori, j, k if at least one of them is active
                // (changed in the previous iteration); otherwise no need to
                // compare i,j,k again.
                if(PV[i].active0 || PV[j].active0 || PV[k].active0 || PV[l].active0){
                  // Set stable = false if swap happens fori,j,k,l
                  if(swap4(i,j,k,l) < 0){ // swap achieved reduction in cost
                    stable4 = false;
                    goto EXIT_i;
                  }
                }
              } //l
            } //k
          } //j
        } //i
      } // while(stable4)
      if(stable4) break; // no swap4 step succeeded; stable3 == stable4 == true;

    EXIT_i: // a swap4 step succeeded; stable4 == false solution is
      // perturbed; now do swap3 again
      stable3 = false;
    }
  } // while(!stable3 || !stable4)
  assert(stable3 && stable4);
  print_solution("Algo: opt34");
  return;
}

// Hill climbing with swap4 (similar to opt4), except a fixed number of
// random swap4 steps are tried.
void opt4_rand(){
  for(size_t l = 0; l < 1000; l++){
    unsigned int i1,i2,i3,i4;
    bool found = false;
    // find 4 distinct random numbers in [0,NP) (uniformly)
    while(!found){
      i1 = random_int(0,NP);
      while((i2 = random_int(0,NP)) == i1);
      while(((i3 = random_int(0,NP)) == i1) || ((i3 = random_int(0,NP)) == i2));
      while(((i4 = random_int(0,NP)) == i1) || ((i4 = random_int(0,NP)) == i2) || ((i4 = random_int(0,NP)) == i3));
      found = true;
    }
    int delta = swap4(i1,i2,i3,i4);
  }
  return;
  print_solution("Algo: opt4_random");
}

// Simulated annealing with swap4.
void opt4_sim_anneal(){
  auto best_cost = total_cost();
  for(size_t k = 1; k < 2000; ++k){ // probabilities
    auto accept_prob = 1.0/((double) k);
    for(size_t l = 0; l < 100; l++){
      // find 3 distinct random numbers in [0,NP) (uniformly)
      unsigned int i1,i2,i3,i4;
      bool found = false;
      while(!found){
        i1 = random_int(0,NP);
        i2 = random_int(0,NP);
        while(i2 == i1){
          i2 = random_int(0,NP);
        }
        i3 = random_int(0,NP);
        while((i3 == i1) || (i3 == i2)){
          i3 = random_int(0,NP);
        }
        i4 = random_int(0,NP);
        while((i4 == i1) || (i4 == i2) || (i4 == i3)){
          i4 = random_int(0,NP);
        }
        found = true;
      }
      int delta = swap4(i1,i2,i3,i4, accept_prob);
    }
    auto tc = total_cost();
    if(tc < best_cost){
      best_cost = tc;
      print_solution("Algo: opt4_sim_anneal");
	  std::cout << "Iteration: " << k << "\n";
      // std::cout << "Iteration: " << k << "\n"
	  // 	<< "Acceptance probability = " << accept_prob << "\n";
    }
  }
  return;
}

// Simulated annealing with swap3.
void opt3_sim_anneal(){
  auto best_cost = total_cost();
  for(size_t k = 1; k < 6000; ++k){ // probabilities
    auto accept_prob = 1.0/((double) k);
    //  double accept_prob = 1.0/sqrt((double)k);
    for(size_t l = 0; l < k+1000; l++){
      // find 3 distinct random numbers in [0,NP) (uniformly)
      unsigned int i1,i2,i3;
      bool found = false;
      while(!found){
        i1 = random_int(0,NP);
        i2 = random_int(0,NP);
        while(i2 == i1){
          i2 = random_int(0,NP);
        }
        i3 = random_int(0,NP);
        while((i3 == i1) || (i3 == i2)){
          i3 = random_int(0,NP);
        }
        found = true;
      }
      int delta = swap3(i1,i2,i3, accept_prob);
    }
    auto tc = total_cost();
    if(tc < best_cost){
      best_cost = tc;
      // print_solution("Algo: opt3_sim_anneal");
      // std::cout << "Iteration: " << k << "\n";
      // std::cout << "Iteration: " << k << "\n"
	  // 	<< "Acceptance probability = " << accept_prob << "\n";
    }
  }
  return;
}


// Hill climbing with swap3 (similar to opt3), except a fixed number of
// random swap3 steps are tried. The loop is broken into two loops so as
// to print the solution in between. This function is almost identical
// to opt3_sim_anneal, except that this one passes 0 probability of
// accepting a bad solution.

void opt3_random(){
  auto best_cost = total_cost();
  for(size_t k = 1; k < 2000; ++k){ // probabilities
    for(size_t l = 0; l < 200; l++){
      // find 3 distinct random numbers in [0,NP) (uniformly)
      unsigned int i1,i2,i3;
      bool found = false;
      while(!found){
        i1 = random_int(0,NP);
        i2 = random_int(0,NP);
        while(i2 == i1){
          i2 = random_int(0,NP);
        }
        i3 = random_int(0,NP);
        while((i3 == i1) || (i3 == i2)){
          i3 = random_int(0,NP);
        }
        found = true;
      }
      int delta = swap3(i1,i2,i3);
    }
    auto tc = total_cost();
    if(tc < best_cost){
      best_cost = tc;
      print_solution("Algo: opt3_random");
      std::cout << "Iteration: " << k << "\n";
    }
  }
  return;
}

// Same as the other version of opt3_rand, but here all iterations are
// done in one loop. Perhaps not very useful if you want to get some
// intermediate results.
void opt3_rand(){
  for(size_t k = 1; k < 800000; ++k){
    // find 3 distinct random numbers in [0,NP) (uniformly)
    unsigned int i1,i2,i3;
    bool found = false;
    while(!found){
      i1 = random_int(0,NP);
      i2 = random_int(0,NP);
      while(i2 == i1){
        i2 = random_int(0,NP);
      }
      i3 = random_int(0,NP);
      while((i3 == i1) || (i3 == i2)){
        i3 = random_int(0,NP);
      }
      found = true;
    }
    int delta = swap3(i1,i2,i3);
    if(delta < 0){
      print_solution("Algo: opt3_random");
      std::cout << "Iteration: " << k << "\n";
    }
  }
  return;
}

// Simulated annealing with swap2.
void opt2_sim_anneal(){
  auto best_cost = total_cost();
  for(size_t k = 1; k < 4000; ++k){ // probabilities
    auto accept_prob = 1.0/((double) k);
    for(size_t l = 0; l < 400; l++){
      // find 3 distinct random numbers in [0,NP) (uniformly)
      unsigned int i1,i2;
      bool found = false;
      while(!found){
        i1 = random_int(0,NP);
        i2 = random_int(0,NP);
        while(i2 == i1){
          i2 = random_int(0,NP);
        }
        found = true;
      }
      int delta = swap2(i1,i2, accept_prob);
    }
    auto tc = total_cost();
    if(tc < best_cost){
      best_cost = tc;
      print_solution("Algo: opt2_sim_anneal");
	  std::cout << "Iteration: " << k << "\n";
      // std::cout << "Iteration: " << k << "\n"
	  // 	<< "Acceptance probability = " << accept_prob << "\n";
    }
  }
  return;
}

// Parse two word options like -f <filename> (stack overflow 865668)
char* get_cmd_option(char ** begin, char ** end, const std::string & option){
  char ** itr = std::find(begin, end, option);
  if (itr != end && ++itr != end)
    {
      return *itr;
    }
  return 0;
}
// Parses single word options like -h for help (stack overflow 865668)
bool cmd_option_exists(char** begin, char** end, const std::string& option)
{
  return std::find(begin, end, option) != end;
}

int main(int argc, char * argv[]){
  // Help menu
  if(cmd_option_exists(argv, argv+argc, "-h")){
    std::cout << "Usage: assign -p <prefs_file> -d <disciplinas_file> -z <seed> OR assign -h (for help)"
              << std::endl
              << "The seed is optional, but useful to get the same solution again"
              << std::endl;
    return 0;
  }
  char *df  = get_cmd_option(argv, argv + argc, "-d"); // disciplinas
  char *pf  = get_cmd_option(argv, argv + argc, "-p"); // profs' preferences
  if(!pf || !df){
    std::cout << "Usage: assign -p <prefs_file> -d <disciplinas_file> -z <seed> OR assign -h (for help)"
              << std::endl
              << "The seed is optional, but useful to get the same solution again"
              << std::endl;
    return 1;
  }

  int seed = 0;
  if(cmd_option_exists(argv, argv+argc, "-z")){
    seed = atoi(get_cmd_option(argv, argv + argc, "-z"));
  }
  randomize(seed);

  std::ifstream fin; // for parsing input files

  // Parse disciplinas
  fin.open(df);
  std::string line;
  uint8_t snd = 0; // short name of a disciplina
  while (std::getline(fin, line)) {
    std::istringstream iss(line);
    Disc d;
    if (!(iss >> d.id >> d.code >> d.nt >> d.nh)) { break; } // error
    d.free = d.nt;
    d.snd = snd;
    DV.push_back(d);
    snd++;
    ND++;
  }
  fin.close();

#if 0
  // test the input disciclinas
  for(auto it = DV.begin(); it != DV.end(); ++it){
    std::cout << (*it).id << ", " << (*it).code << ", " << (*it).nt << ", " << (*it).nh << "\n";
  }
#endif

  // Make a map of <disciplina, index>
  for(auto it = DV.begin() ; it != DV.end(); ++it){
    // Make string of id+code, e.g., id=20, code=F will give string "20F"
    std::ostringstream ost;
    ost << (*it).id << (*it).code;
    std::string s = ost.str();
    // Put <s,index of it> in the map
    DM[s] = it-DV.begin();
  }

  // Parse preferences of profs
  fin.open(pf);
  std::istringstream iss1; // for lines
  std::istringstream iss2; // for tokens in lines
  NP = 0;
  uint8_t snp = 0; // short name of a prof
  while (std::getline(fin, line)) {
    iss1.str(line);
    Prof pi;
    // tokenise iss1
    std::string token;
    // first token is name
    std::getline(iss1,token, ',');
    iss2.str(token);
    if(!(iss2 >> pi.name)) break;
    iss2.clear();
    // second token is min. number of hours teaching
    std::getline(iss1,token, ',');
    iss2.str(token);
    if(!(iss2 >> pi.hmin)) break;
    iss2.clear();
    // third token is max. number of hours teaching
    std::getline(iss1,token, ',');
    iss2.str(token);
    if(!(iss2 >> pi.hmax)) break;
    iss2.clear();
    pi.snp = snp;
    // Remaining tokens are preferences; each preference is 3 words.
    unsigned int num_prefs = 0;
    while(std::getline(iss1,token, ',')){
      iss2.str(token);
      unsigned int s0;
      std::string s1;
      std::string s2;
      std::string s3;
      iss2 >> s0 >> s1 >> s2 >> s3;
      std::array<uint8_t,4> x = {DM[s1],DM[s2],DM[s3],pi.snp};
      Assignment ass(x);
	  pi.prefs.push_back(ass);
      cost_map[ass] = s0;
      iss2.clear();
      num_prefs++;
    }
    iss1.clear();
    if(num_prefs == 0) pi.min_pref_cost = 0;
    PV.push_back(pi);
    NP++;
    snp++;
  }
  fin.close();

  preprocess();

  // Possible scripts to use the algorithms
#if 1
  // Simulated annealing with opt3_sim_anneal
  init();
  // print_solution("Initial solution");
  opt3_sim_anneal();
  print_solution("Algo: opt3_sim_anneal");
  // opt2_sim_anneal();
  // opt4();
  // print_solution("Algo: opt4");
#endif

#if 0
  init();
  opt4();
  print_solution("Algo: opt4");
#endif

#if 0
  // Try opt2. If the cost is lower than a threshold, then try opt3; if
  // the cost after opt3 is lower than a threshold, then try opt4.
  // reasonable value of the objective function, and then use opt3 on it
  unsigned int i = 10000;
  unsigned int tc2 = 600;
  unsigned int tc3 = 130;
  while(i-- > 0){
    init();
    opt2();
    if(total_cost() > tc2) continue;
    opt3();
    if(total_cost() > tc3) continue;
    print_solution("Algo: opt3");
    opt4();
    print_solution("Algo: opt4");
  }
#endif
  return 0;
}
