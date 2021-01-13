//  This script is adapted from Bellman's GAP example scripts
//  http://gapc.eu/
//  Copyright (C) Stefan Janssen, Georg Sauthoff, Robert Giegerich.
//    Technische Fakultät, Bielefeld University, Germany.
//  Files:
//    hishapes_ambiguous.gap
//    hishapes_unambiguous.gap
//  License: GPL
  

import rna
import filter_add

import mfe_filter
import hishape_mfe_filter
import mfe_answer
import mfe_answer_v4
import pf_answer
import pf_filter
import hishape_answer_v4
import hishape_filter
input rna

type pfanswer = extern
type mfeanswer = extern
type mfeanswer_v4 = extern
type Rope = extern
type hishapeanswer_v4 = extern
type base_t = extern

signature Canonical_Algebra(alphabet,answer) {
	answer sadd(Subsequence,answer);
	answer cadd(answer,answer);
	answer cadd_Pr(answer,answer);
	answer cadd_Pr_Pr(answer,answer);
	answer cadd_Pr_Pr_Pr(answer,answer);
	answer ambd(answer,Subsequence,answer);
	answer ambd_Pr(answer,Subsequence,answer);
	answer nil(Subsequence);
	answer nil_Pr(Subsequence);
	answer edl(Subsequence,answer);
	answer edr(answer,Subsequence);
	answer edlr(Subsequence,answer,Subsequence);
	answer drem(answer);
	answer is(answer);
	answer sr(Subsequence,answer,Subsequence);
	answer hl(Subsequence,Subsequence,Subsequence);
	answer sp(Subsequence,answer,Subsequence);
	answer bl(Subsequence,answer);
	answer br(answer,Subsequence);
	answer il(Subsequence,answer,Subsequence);
	answer ml(Subsequence,answer,Subsequence);
	answer mldr(Subsequence,answer,Subsequence,Subsequence);
	answer mladr(Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldlr(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mladlr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldladr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladldr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldl(Subsequence,Subsequence,answer,Subsequence);
	answer mladl(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer addss(answer,Subsequence);
	answer ssadd(Subsequence,answer);
	answer trafo(answer);
	answer incl(answer);
	answer combine(answer,answer);
	answer acomb(answer,Subsequence,answer);
	choice [answer] h([answer]);
}

algebra mfeV4 implements Canonical_Algebra(alphabet = char, answer = mfeanswer_v4) {
	mfeanswer_v4 sadd(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.subword.i = lb.i;

		return res;
	}

	mfeanswer_v4 cadd(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

                res.outermost_initstem_i = le.outermost_initstem_i;
                res.next_ss_i = le.outermost_initstem_j;
                if (re.outermost_initstem_j == -1) {  // if the right part is empty, happen when sadd(nil) in struct 
                  res.outermost_initstem_j = le.outermost_initstem_j;
                  res.next_ss_j = re.firstStem.j;
                } else {
                  res.outermost_initstem_j = re.outermost_initstem_j;
                  res.next_ss_j = re.outermost_initstem_i;
                }
               
		return res;
	}

	mfeanswer_v4 cadd_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 ambd(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 ambd_Pr(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 nil(Subsequence loc) {
		mfeanswer_v4 res;          
		
		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;

                res.outermost_initstem_i = -1;
                res.outermost_initstem_j = -1;
                res.next_ss_i = -1;
                res.next_ss_j = -1;

		return res;
	}

	mfeanswer_v4 nil_Pr(Subsequence loc) {
		mfeanswer_v4 res;           
		
		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;

		return res;
	}

	mfeanswer_v4 edl(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;

		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
                       
		return res;
	}

	mfeanswer_v4 edr(mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;

		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem);
		res.subword.j = rb.j;

		return res;
	}

	mfeanswer_v4 edlr(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;

		res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem);
		res.subword.i = lb.i;
		res.subword.j = rb.j;

		return res;
	}

	mfeanswer_v4 drem(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;

		return res;
	}

	mfeanswer_v4 is(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + termau_energy(e.firstStem, e.firstStem);

                res.outermost_initstem_i = e.firstStem.i;
                res.outermost_initstem_j = e.firstStem.j;
                res.next_ss_i = e.firstStem.j;
                res.next_ss_j = e.firstStem.j;

		return res;
	}

	mfeanswer_v4 sr(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		

		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = false;
		
		return res;
	}

	mfeanswer_v4 hl(Subsequence lb,Subsequence region,Subsequence rb) {
		mfeanswer_v4 res;                // based on empty
		          
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = hl_energy(region);  // Note: using the positions before closing bp
		res.lastStem = res.firstStem;    // lastStem and subword have the same content of firstStem
		res.subword = res.firstStem;     // lastStem and subword have the same content of firstStem

                res.prevE = 0;
                res.prevIsMl = false;
		
		return res;
	}

	mfeanswer_v4 sp(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = e.energy;  // this is correct, bl_energy in bl, sp without energy
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

		return res;
	}

	mfeanswer_v4 bl(Subsequence lregion,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		Subsequence innerstem;
		innerstem.seq = lregion.seq;
		innerstem.i = lregion.i-1;
		innerstem.j = e.firstStem.j+1;
		// bl_energy(..) includes the termau_energy ?! Question 1
		res.energy = e.energy + bl_energy(lregion,innerstem);  // use lregion and innerstem with closing bp
		res.subword.i = lregion.i;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = false; 
		
		return res;
	}

	mfeanswer_v4 br(mfeanswer_v4 e,Subsequence rregion) {
		mfeanswer_v4 res = e;
		              
		Subsequence innerstem;
		innerstem.seq = rregion.seq;
		innerstem.i = e.firstStem.i-1;
		innerstem.j = rregion.j+1;
		res.energy = e.energy + br_energy(innerstem, rregion);  
		res.subword.j = rregion.j;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = false; 

		return res;
	}

	mfeanswer_v4 il(Subsequence lregion,mfeanswer_v4 e,Subsequence rregion) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + il_energy(lregion, rregion);
		res.subword.i = lregion.i;
		res.subword.j = rregion.j;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = false; 

		return res;
	}

	mfeanswer_v4 ml(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;

		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
                // ml_energy() for penalty of creating a multiloop
                // ul_energy() is for how many hairpin loops within multiloop
		res.energy = ml_energy() + ul_energy() + e.energy + termau_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

                // (...(...(...)..(...)..)...) since both addss() and sadd() do not change the energies, for late 
                res.prevSubword = e.subword;  // i and j positions of ...(...)..(...).., not (...)..(...), but does not matter since we do not need the infos
                res.prevE = e.energy;  // energy of ...(...)..(...).. == energy of (...)..(...)
                res.prevIsMl = true; 
		
		return res;
	}

	mfeanswer_v4 mldr(Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = true; 
		
		return res;
	}

	mfeanswer_v4 mladr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldlr(Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = true; 
		
		return res;
	}

	mfeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldl(Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(res.firstStem,res.firstStem) + termau_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

                res.prevSubword = e.subword;
                res.prevE = e.energy;
                res.prevIsMl = true; 
		
		return res;
	}

	mfeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;

		return res;
	}

	mfeanswer_v4 addss(mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + ss_energy(rb);
		res.subword.j = rb.j;
  
                res.outermost_initstem_i = e.firstStem.i;
                res.outermost_initstem_j = e.firstStem.j;
                res.next_ss_i = rb.i;            // addss() in ml_comps
                res.next_ss_j = rb.j;     // is() initilize the value with rightmost values (...)[i,j]
                                          // addss() extends the next_ss with (...)[...]
                                          // cadd() update values 
                                          // nil() marks the empty-tags

		return res;
	}

	mfeanswer_v4 ssadd(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v4 trafo(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		return res;
	}

	mfeanswer_v4 incl(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = ul_energy() + e.energy;
		
		return res;
	}

	mfeanswer_v4 combine(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v4 acomb(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	choice [mfeanswer_v4] h([mfeanswer_v4] i) {
		return list(minimum(i));
	}
}






algebra count auto count ;

algebra enum auto enum ;


algebra pretty implements Canonical_Algebra(alphabet = char, answer = string) {
	string sadd(Subsequence lb,string e) {
		string res;
		append(res, '.');
		append(res, e);
		return res;
	}

	string cadd(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string cadd_Pr_Pr_Pr(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string ambd(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, '.');
		append(res, re);
		return res;
	}

	string ambd_Pr(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, '.');
		append(res, re);
		return res;
	}

	string nil(Subsequence loc) {
		string r;
		return r;
	}

	string nil_Pr(Subsequence loc) {
		string r;
		return r;
	}

	string edl(Subsequence lb,string e) {
		string res;
		append(res, '.');
		append(res, e);
		return res;
	}

	string edr(string e,Subsequence rb) {
		string res;
		append(res, e);
		append(res, '.');
		return res;
	}

	string edlr(Subsequence lb,string e,Subsequence rb) {
		string res;
		append(res, '.');
		append(res, e);
		append(res, '.');
		return res;
	}

	string drem(string e) {
		return e;
	}

	string is(string e) {
		return e;
	}

	string sr(Subsequence lb,string e,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, e);
		append(res, ')');
		return res;
	}

	string hl(Subsequence lb,Subsequence region,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, '.', size(region));
		append(res, ')');
		return res;
	}

	string sp(Subsequence lb,string e,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, e);
		append(res, ')');
		return res;
	}

	string bl(Subsequence lregion,string e) {
		string res;
		append(res, '.', size(lregion));
		append(res, e);
		return res;
	}

	string br(string e,Subsequence rregion) {
		string res;
		append(res, e);
		append(res, '.', size(rregion));
		return res;
	}

	string il(Subsequence lregion,string e,Subsequence rregion) {
		string res;
		append(res, '.', size(lregion));
		append(res, e);
		append(res, '.', size(rregion));
		return res;
	}

	string ml(Subsequence lb,string e,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, e);
		append(res, ')');
		return res;
	}

	string mldr(Subsequence lb,string e,Subsequence dr,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, e);
		append(res, '.');
		append(res, ')');
		return res;
	}

	string mladr(Subsequence llb,Subsequence lb,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string mldlr(Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, ')');
		return res;
	}

	string mladlr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string mldladr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string mladldr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
		return res;
	}

	string mldl(Subsequence lb,Subsequence dl,string e,Subsequence rb) {
		string res;
		append(res, '(');
		append(res, '.');
		append(res, e);
		append(res, ')');
		return res;
	}

	string mladl(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, "))", 2);
		return res;
	}

	string addss(string e,Subsequence rb) {
		string res;
		append(res, e);
		append(res, '.', size(rb));
		return res;
	}

	string ssadd(Subsequence lb,string e) {
		string res;
		append(res, '.', size(lb));
		append(res, e);
		return res;
	}

	string trafo(string e) {
		return e;
	}

	string incl(string e) {
		return e;
	}

	string combine(string le,string re) {
		string res;
		append(res, le);
		append(res, re);
		return res;
	}

	string acomb(string le,Subsequence b,string re) {
		string res;
		append(res, le);
		append(res, '.');
		append(res, re);
		return res;
	}

	choice [string] h([string] i) {
		//~ return list(minimum(i));
		//return i;
                return unique(i);
	}
}



///////////////////////////////////////////////////////////////////////////
///////////////////////////// hishapeb ////////////////////////////////////
algebra hishapeh implements Canonical_Algebra(alphabet = char, answer = Rope) {
	Rope sadd(Subsequence b,Rope e) {
		Rope emptyShape;
		Rope res;
		
		if (e == emptyShape) {
			append(res, '_');
			append(res, e);
			return res;
		} else {
			return e;
		}
	}

	Rope cadd(Rope le,Rope re) {	
		if (re == "_") {
			return le;
		} else {
			return le + re;
		}
	}

	Rope cadd_Pr(Rope le,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope cadd_Pr_Pr(Rope le,Rope re) {
		if (re == "_") {
			return le;
		} else {
			return le + re;
		}
	}

	Rope cadd_Pr_Pr_Pr(Rope le,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope ambd(Rope le,Subsequence b,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope ambd_Pr(Rope le,Subsequence b,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope nil(Subsequence loc) {
		Rope r;
		return r;
	}

	Rope nil_Pr(Subsequence loc) {
		Rope r;
		return r;
	}

	Rope edl(Subsequence lb,Rope e) {
		return e;
	}

	Rope edr(Rope e,Subsequence rb) {
		return e;
	}

	Rope edlr(Subsequence lb,Rope e,Subsequence rb) {
		return e;
	}

	Rope drem(Rope e) {
		return e;
	}

	Rope is(Rope e) {
		return e;
	}

	Rope sr(Subsequence lb,Rope e,Subsequence rb) {
		return e;
	}

	Rope hl(Subsequence lb,Subsequence region,Subsequence rb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
		append(res, ',');
		return res;
	}

	Rope sp(Subsequence lb,Rope e,Subsequence rb) {
		return e;
	}

	Rope bl(Subsequence lregion,Rope e) {
		return e;
	}

	Rope br(Rope e,Subsequence rregion) {
		return e;
	}

	Rope il(Subsequence lregion,Rope e,Subsequence rregion) {
		return e;
	}

	Rope ml(Subsequence lb,Rope e,Subsequence rb) {
		return e;
	}

	Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
		return e;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
		return e;
	}

	Rope mladlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mldladr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mladldr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
		return e;
	}

	Rope mladl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope addss(Rope e,Subsequence rb) {
		return e;
	}

	Rope ssadd(Subsequence lb,Rope e) {
		return e;
	}

	Rope trafo(Rope e) {
		return e;
	}

	Rope incl(Rope e) {
		return e;
	}

	Rope combine(Rope le,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	Rope acomb(Rope le,Subsequence b,Rope re) {
		Rope res;
		append(res, le);
		append(res, re);
		return res;
	}

	choice [Rope] h([Rope] i) {
		return unique(i);
	}
}

algebra hishapehplus extends hishapeh {
	Rope ml(Subsequence lb,Rope e,Subsequence rb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mladlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mldladr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mladldr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mladl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}
}


algebra hishapem extends hishapeh {
	Rope ml(Subsequence lb,Rope e,Subsequence rb) {
		Rope res;       
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mldr(Subsequence lb,Rope e,Subsequence dr,Subsequence rb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mldlr(Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mladlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mldladr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mladldr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mldl(Subsequence lb,Subsequence dl,Rope e,Subsequence rb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}

	Rope mladl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
                append(res, "m,(,", 4);
                append(res, e);
                append(res, "),", 2);
		return res;
	}
}


algebra hishapeb extends hishapem {
	Rope sp(Subsequence lb,Rope e,Subsequence rb) {
		Rope res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) append(res, ".5", 2);
		append(res, e);
		return res;
	}

	Rope bl(Subsequence lregion,Rope e) {
                Rope res;
		append(res, "b,", 2);
		append(res, e);
		return res;
	}

	Rope br(Rope e,Subsequence rregion) {
                Rope res;
		append(res, "b,", 2);
		append(res, e);
		return res;
	}

	Rope il(Subsequence lregion,Rope e,Subsequence rregion) {
                Rope res;
		append(res, "i,", 2);
		append(res, e);
		return res;
	}
}


grammar gra_nodangle uses Canonical_Algebra(axiom = struct) {
  tabulated {
    hairpin, iloop, closed, ml_comps, ml_comps1, struct }  // leftB, rightB,  stack, 

  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle = drem(initstem) # h;
          
  initstem = is(closed) # h;  // for adding termau_energy last position vs. incl() for energy initiating stem ul_energy()

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} # h;

  stack     = sr(BASE,                                closed,                                BASE) with stackpairing # h;
  hairpin   = hl(BASE,                          REGION with minsize(3),          BASE) with basepairing # h;
  leftB     = sp(BASE, bl(REGION, closed), BASE)                                            with basepairing # h;
  rightB    = sp(BASE, br(closed, REGION), BASE)                                            with basepairing # h;
  iloop = sp(BASE, il(REGION with maxsize(30), closed, REGION with maxsize(30)), BASE) with basepairing # h;

  multiloop = ml(BASE,                          ml_comps,                        BASE) with basepairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;  // incl is for ul_energy()

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}


grammar gra_dangle uses Canonical_Algebra(axiom = struct) {
  tabulated {
    hairpin, iloop, closed, ml_comps, ml_comps1, struct }  // leftB, rightB,  stack, 

  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = edl (BASE, initstem) |
              edr (initstem, BASE) |             
              edlr(BASE, initstem, BASE) |        
              drem(initstem) # h;                
  initstem = is(closed) # h;

  closed    = {stack                        | 
               hairpin                      |
               leftB                        | 
               rightB                       | 
               iloop                        | 
               multiloop} # h;

  stack     = sr(BASE,                                closed,                                BASE) with stackpairing # h;
  hairpin   = hl(BASE,                          REGION with minsize(3),          BASE) with basepairing # h;
  leftB     = sp(BASE, bl(REGION, closed), BASE)                                            with basepairing # h;
  rightB    = sp(BASE, br(closed, REGION), BASE)                                            with basepairing # h;
  iloop = sp(BASE, il(REGION with maxsize(30), closed, REGION with maxsize(30)), BASE) with basepairing # h;

  multiloop = ml(BASE,                          ml_comps,                        BASE) with basepairing |
              mldl(BASE, BASE,                  ml_comps,                        BASE) with basepairing |
              mldr(BASE,                        ml_comps,                  BASE, BASE) with basepairing |
              mldlr(BASE, BASE,                 ml_comps,                  BASE, BASE) with basepairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}


// gapc ~.gap -i pp_mfeV4 -o pp_mfeV4_.cc --kbest *
instance pp_mfeV4 = gra_nodangle (pretty * mfeV4);
// gapc ~.gap -i pp_mfeV4_d1 -o pp_mfeV4_d1_.cc --kbest *
instance pp_mfeV4_d1 = gra_dangle (pretty * mfeV4);

// gapc ~.gap -i rope_mfeV4_d1 -o rope_mfeV4_d1_.cc --kbest *
instance rope_mfeV4_d1 = gra_nodangle (hishapeh * mfeV4);