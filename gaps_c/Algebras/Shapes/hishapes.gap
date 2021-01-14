import rna
import filter_add
//##import algebra_add
//import list_add

import mfe_filter
import hishape_mfe_filter
import mfe_answer
////import mfe_answer_v2
import mfe_answer_v4
import pf_answer
import pf_filter
////import hishape_answer_v2
////import hishape_answer_v3
import hishape_answer_v4
import hishape_filter
input rna


type pfanswer = extern
type mfeanswer = extern  // change in accordance with 'import mfe_answer'
////////type mfeanswer = (int energy, Subsequence firstStem, Subsequence lastStem)
//##type mfeanswer = (int energy)
////////type mfeanswer_dbg = (int energy, Subsequence firstStem, Subsequence lastStem, string rep)
////type mfeanswer_v2 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
////type mfeanswer_v2 = extern
////type mfeanswer_v3 = (int energy, Subsequence firstStem, Subsequence lastStem, string rep, string pp)
//////type mfeanswer_v4 = (int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep, string pp)
type mfeanswer_v4 = extern

type Rope = extern
//type hishapeanswer_v2 = (Rope rope, int energy, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
////type hishapeanswer_v2 = extern
//type hishapeanswer_v3 = (Rope rope, Subsequence firstStem, Subsequence lastStem, Subsequence subword, string rep)
////type hishapeanswer_v3 = extern
type hishapeanswer_v4 = extern

type base_t = extern

// answer files can be extracted from generated header files directly



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
	answer hl(Subsequence,Subsequence,Subsequence,Subsequence,Subsequence);
	answer sp(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer bl(Subsequence,answer);
	answer br(answer,Subsequence);
	answer il(Subsequence,answer,Subsequence);
	answer ml(Subsequence,Subsequence,answer,Subsequence,Subsequence);
	answer mldr(Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladr(Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldlr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladlr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldladr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mladldr(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence,Subsequence);
	answer mldl(Subsequence,Subsequence,Subsequence,answer,Subsequence,Subsequence);
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
		
		Rope empty_shape;
		if (e.hi_rope == empty_shape ) {
			res.hi_rope = "_";
                        res.hi_i_j_rope = "_";
		}
                res.before_is = false;

		res.subword.i = lb.i;


		string o6;
		append(o6, '.');
		append(o6, e.pp);
		res.pp = o6;

		string o;
		append(o, "sadd{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 cadd(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		if (re.hi_rope == "_" ) {
			res.hi_rope = le.hi_rope;
                        res.hi_i_j_rope = le.hi_i_j_rope;
		} else {
			res.hi_rope = le.hi_rope + re.hi_rope;
                        res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
		}
                res.before_is = false;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		string o;
		append(o, "cadd{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, le.rep);
		append(o, ",", 1);
		append(o, re.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 cadd_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_rope = le.hi_rope + re.hi_rope;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
                res.before_is = false;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		string o;
		append(o, "cadd'{", 6);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, le.rep);
		append(o, ",", 1);
		append(o, re.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		if (re.hi_rope == "_" ) {
			res.hi_rope = le.hi_rope;
                        res.hi_i_j_rope = le.hi_i_j_rope;
		} else {
			res.hi_rope = le.hi_rope + re.hi_rope;
                        res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
		}
                res.before_is = false;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		string o;
		append(o, "cadd''{", 7);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, le.rep);
		append(o, ",", 1);
		append(o, re.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_rope = le.hi_rope + re.hi_rope;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
                res.before_is = false;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		string o;
		append(o, "cadd'''{", 8);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, le.rep);
		append(o, ",", 1);
		append(o, re.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 ambd(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_rope = le.hi_rope + re.hi_rope;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
                res.before_is = false;

		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;

		string o1;
		string o2;
		append(o1, "ambd{", 5);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, le.rep);
		append(o1, ",", 1);
		append(o1, re.rep);
		append(o1, ",min(dr_energy(", 15);
		append(o2, dr_energy(le.firstStem, le.firstStem));
		append(o2, "),dl_energy(", 12);
		append(o2, ")=",2);
		append(o2, min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)));
		append(o2, ")}", 2);
		string o;
		append(o,o1);
		append(o,o2);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;


		return res;
	}

	mfeanswer_v4 ambd_Pr(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_rope = le.hi_rope + re.hi_rope;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;
                res.before_is = false;

		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;

		string o1;
		string o2;
		append(o1, "ambd'{", 5);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, le.rep);
		append(o1, ",", 1);
		append(o1, re.rep);
		append(o1, ",min(dr_energy(", 15);
		append(o2, dr_energy(le.firstStem, le.firstStem));
		append(o2, "),dl_energy(", 12);
		append(o2, ")=",2);
		append(o2, min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem)));
		append(o2, ")}", 2);
		string o;
		append(o,o1);
		append(o,o2);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 nil(Subsequence loc) {
		mfeanswer_v4 res;          // based on empty
		
		Rope r;                     // r is empty
		res.hi_rope = r;
                res.hi_i_j_rope = r;
                res.before_is = false;

		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;

		string o6;
		res.pp = o6;

		string o;
		append(o, "nil{", 4);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, "0}", 2);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 nil_Pr(Subsequence loc) {
		mfeanswer_v4 res;           // based on empty
		
		Rope r;
		res.hi_rope = r;
                res.hi_i_j_rope = r;
                res.before_is = false;

		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;
		
		string o6;
		res.pp = o6;

		string o;
		append(o, "nil'{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, "0}", 2);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 edl(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
                res.before_is = false;

		Subsequence innerstem;
		innerstem.seq = e.firstStem.seq;
		innerstem.i = e.firstStem.i+1;
		innerstem.j = e.firstStem.j-1;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem);
                res.shouldPush = true;
                if (res.is2BpHL) {
//                  if ((hl_energy(innerstem)+termau_energy(innerstem,innerstem)+ext_mismatch_energy(innerstem,innerstem))<res.energy) {
//                    res.shouldPush = false;
//                  }
//                  if ((hl_energy(e.firstStem)+termau_energy(e.firstStem, e.firstStem)+dl_energy(e.firstStem, e.firstStem))<res.energy) {
//                    res.shouldPush = false;
//                  }
                    if (0<res.energy) {
                      res.shouldPush = false;
                    }
                }   
		res.subword.i = lb.i;
             
                

		string o6;
		append(o6, '.');
		append(o6, e.pp);
		res.pp = o6;

		string o;
		append(o, "edl{", 4);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, dl_energy(e.firstStem, e.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 edr(mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		res.before_is = false;

		Subsequence innerstem;
		innerstem.seq = e.firstStem.seq;
		innerstem.i = e.firstStem.i+1;
		innerstem.j = e.firstStem.j-1;
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem);
                res.shouldPush = true;
                if (res.is2BpHL) {
//                  if ((hl_energy(innerstem)+termau_energy(innerstem,innerstem)+ext_mismatch_energy(innerstem,innerstem))<res.energy) {
//                    res.shouldPush = false;
//                  }
//                  if ((hl_energy(e.firstStem)+termau_energy(e.firstStem, e.firstStem)+dr_energy(e.firstStem, e.firstStem))<res.energy) {
//                    res.shouldPush = false;
//                  }
                    if (0<res.energy) {
                      res.shouldPush = false;
                    }
                }   
		res.subword.j = rb.j;

		string o6;
		append(o6, e.pp);
		append(o6, '.');
		res.pp = o6;

		string o;
		append(o, "edr{", 4);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, dr_energy(e.firstStem, e.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 edlr(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		res.before_is = false;

		Subsequence innerstem;
		innerstem.seq = e.firstStem.seq;
		innerstem.i = e.firstStem.i+1;
		innerstem.j = e.firstStem.j-1;
		res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem);
                res.shouldPush = true;
                if (res.is2BpHL) {
//                  if ((hl_energy(innerstem)+termau_energy(innerstem,innerstem)+ext_mismatch_energy(innerstem,innerstem))<res.energy) {
//                    res.shouldPush = false;
//                  }
//                  if ((hl_energy(e.firstStem)+termau_energy(e.firstStem, e.firstStem)+ext_mismatch_energy(e.firstStem, e.firstStem))<res.energy) {
//                    res.shouldPush = false;
//                  }
                    if (0<res.energy) {
                      res.shouldPush = false;
                    }
                }   
		res.subword.i = lb.i;
		res.subword.j = rb.j;

		string o6;
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		res.pp = o6;

		string o;
		append(o, "edlr{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ext_mismatch_energy(e.firstStem, e.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 drem(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
                res.before_is = false;

		Subsequence innerstem;
		innerstem.seq = e.firstStem.seq;
		innerstem.i = e.firstStem.i+1;
		innerstem.j = e.firstStem.j-1;
                //DO NOT ENERGY CALCULATION    res.energy = e.energy;
                res.shouldPush = true;
                if (res.is2BpHL) {
//                  if ((hl_energy(innerstem)+termau_energy(innerstem,innerstem)+ext_mismatch_energy(innerstem,innerstem))<res.energy) {
//                    res.shouldPush = false;
//                  }
//                  if ((hl_energy(e.firstStem)+termau_energy(e.firstStem, e.firstStem))<res.energy) {
//                    res.shouldPush = false;
//                  }
                    if (0<res.energy) {
                      res.shouldPush = false;
                    }
                } 
                //DO NOT subword REDEFINITION

		res.pp = e.pp;

		string o;
		append(o, "drem{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, e.pp);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 is(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		Rope hi_i_j_res;
                append(hi_i_j_res, e.hi_i_j_rope);

                append(hi_i_j_res, '_');
                append(hi_i_j_res, e.i1);
                append(hi_i_j_res, '_');
                append(hi_i_j_res, e.j1);
                append(hi_i_j_res, ',');
                res.hi_i_j_rope = hi_i_j_res;
                res.before_is = false;

		res.energy = e.energy + termau_energy(e.firstStem, e.firstStem);
		
		res.pp = e.pp;

		string o;
		append(o, "is{", 3);
		append(o, e.rep);
		append(o, ",termau_energy(",15);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, termau_energy(e.firstStem, e.firstStem));
		append(o, ")}", 2);
                append(o, "[", 1);
                append(o, e.pp);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 sr(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
                res.i1 = lb.i;
                res.j1 = rb.j;
                res.is2BpHL = false;
                res.is2BpML = false;

		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, '(');
		append(o6, e.pp);
		append(o6, ')');
		res.pp = o6;

		string o;
		append(o, "sr{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, sr_energy(res.firstStem,res.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res;                // based on empty
		
		//Rope hi_res;
                Rope hi_i_j_res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		//append(hi_res, pos);
                append(hi_i_j_res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) 
                {
                  //append(hi_res, ".5", 2);
                  append(hi_i_j_res, ".5", 2);
                }
		//append(hi_res, ',');
		//res.hi_rope = hi_res;
                res.hi_i_j_rope = hi_i_j_res;

                res.i0 = llb.i;
                res.j0 = rrb.j;
                res.i1 = llb.i;
                res.j1 = rrb.j;
                res.is2BpHL = true;
                res.is2BpML = false;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = hl_energy(region) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;    // lastStem and subword have the same content of firstStem
		res.subword = res.firstStem;     // lastStem and subword have the same content of firstStem
		
		string o6;
		append(o6, "((",2);
		append(o6, '.', size(region));
		append(o6, "))",2);
		res.pp = o6;

		string o;
		append(o, "hl{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, hl_energy(region) + sr_energy(res.firstStem,res.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 sp(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.i1 = llb.i;
                res.j1 = rrb.j;
                res.is2BpHL = false;
                res.is2BpML = false;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, "((",2);
		append(o6, e.pp);
		append(o6, "))",2);
		res.pp = o6;

		string o;
		append(o, "sp{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, sr_energy(res.firstStem,res.firstStem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;


		return res;
	}

	mfeanswer_v4 bl(Subsequence lregion,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		Subsequence innerstem;
		innerstem.seq = lregion.seq;
		innerstem.i = lregion.i-1;
		innerstem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(lregion,innerstem);
		res.subword.i = lregion.i;
		
		string o6;
		append(o6, '.', size(lregion));
		append(o6, e.pp);
		res.pp = o6;

		string o;
		append(o, "bl{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, lregion);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, bl_energy(lregion,innerstem));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



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
		
		string o6;
		append(o6, e.pp);
		append(o6, '.', size(rregion));
		res.pp = o6;

		string o;
		append(o, "br{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, br_energy(innerstem, rregion));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 il(Subsequence lregion,mfeanswer_v4 e,Subsequence rregion) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + il_energy(lregion, rregion);
		res.subword.i = lregion.i;
		res.subword.j = rregion.j;
		
		string o6;
		append(o6, '.', size(lregion));
		append(o6, e.pp);
		append(o6, '.', size(rregion));
		res.pp = o6;

		string o;
		append(o, "il{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, il_energy(lregion, rregion));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 ml(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
// ml_energy() for penalty of creating a multiloop
// ul_energy() is for how many hairpin loops within multiloop
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "ml{", 3);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mldr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "mldr{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mladr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

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
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o1;
		string o2;
		string o3;
		append(o1, "mladr{380,", 10);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, e.rep);
		append(o1, ",", 1);
		append(o1, "min(dri_energy(", 15);
		append(o1, dri_energy(innerstem,innerstem));
		append(o2, "),dr_energy(", 12);
		append(o2, dr_energy(e.lastStem, e.lastStem));
		append(o2, ")=",2);
		append(o2, min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)));
		append(o3, "),sr_energy(",12);
		append(o3, sr_energy(res.firstStem,res.firstStem));
		append(o3, "), termau_energy(", 16);
		append(o3, termau_energy(innerstem,innerstem));
		append(o3, ")}", 2);
		string o;
		append(o,o1);
		append(o,o2);
		append(o,o3);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "mldlr{", 6);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

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
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o1;
		string o2;
		string o3;
		string o4;
		string o5;
		append(o1, "mladlr{380,", 11);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, e.rep);
		append(o1, ",min(dli_energy(",16);
		append(o1, dli_energy(innerstem,innerstem));
		append(o2, "),dl_energy(",12);
		append(o2, dl_energy(e.firstStem, e.firstStem));
		append(o2, "=",1);
		append(o2, min(dli_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)));
		append(o3, "),min(dri_energy(",17);
		append(o3, dri_energy(innerstem,innerstem));
		append(o3, "),dl_energy(",12);
		append(o3, dr_energy(e.lastStem, e.lastStem));
		append(o4, "=",1);
		append(o4, min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)));
		append(o4, "),sr_energy(",12);
		append(o4, sr_energy(res.firstStem,res.firstStem));
		append(o5, "),termau_energy(",16);
		append(o5, termau_energy(innerstem,innerstem));
		append(o5, ")}",2);
		string o;
		append(o,o1);
		append(o,o2);
		append(o,o3);
		append(o,o4);
		append(o,o5);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

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
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o1;
		string o2;
		string o3;
		string o4;
		append(o1, "mldladr{380, ", 13);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, e.rep);
		append(o1, "), dli_energy(", 14);
		append(o2, dli_energy(innerstem,innerstem));
		append(o2, "), min(dri_energy(", 18);
		append(o2, dri_energy(innerstem,innerstem));
		append(o2, "), dr_energy(", 13);
		append(o3, dr_energy(e.lastStem,e.lastStem));
		append(o3, "))=", 3);
		append(o3, min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem,e.lastStem)));
		append(o3, ", sr_energy(", 12);
		append(o4, sr_energy(res.firstStem,res.firstStem));
		append(o4, "), termau_energy(", 17);
		append(o4, termau_energy(innerstem,innerstem));
		append(o4, ")}", 2);
		string o;
		append(o, o1);
		append(o, o2);
		append(o, o3);
		append(o, o4);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
		}

	mfeanswer_v4 mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

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
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "mladldr{", 8);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "mldl{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.is2BpHL = false;
                res.is2BpML = true;

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

		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;

		string o;
		append(o, "mladl{", 6);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ml_energy() + ul_energy()+min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem));
		append(o, "}", 1);
                //append(o, "[", 1);
                //append(o, o6);
                //append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 addss(mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + ss_energy(rb);
		res.subword.j = rb.j;

		string o6;
		append(o6, e.pp);
		append(o6, '.', size(rb));
		res.pp = o6;

		string o;
		append(o, "addss{", 6);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ss_energy(rb));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 ssadd(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		res.subword.i = lb.i;
		
		string o6;
		append(o6, '.', size(lb));
		append(o6, e.pp);
		res.pp = o6;

		string o;
		append(o, "ssadd{40,", 9);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, ss_energy(lb));
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 trafo(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.pp = e.pp;

		string o;
		append(o, "trafo{", 6);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, e.pp);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 incl(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = ul_energy() + e.energy;
		
		res.pp = e.pp;

		string o;
		append(o, "incl{", 5);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, e.rep);
		append(o, ",", 1);
		append(o, "40}", 3);
                append(o, "[", 1);
                append(o, e.pp);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 combine(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		Rope lr = le.hi_rope;
		Rope rr = re.hi_rope;
		Rope hi_res;
		append(hi_res, lr);
		append(hi_res, rr);
		res.hi_rope = hi_res;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		string o;
		append(o, "combine{", 8);
		//append(o, res.firstStem);
		//append(o, res.subword);
		//append(o, res.lastStem);
		append(o, le.rep);
		append(o, ",", 1);
		append(o, re.rep);
		append(o, "}", 1);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;



		return res;
	}

	mfeanswer_v4 acomb(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		Rope lr = le.hi_rope;
		Rope rr = re.hi_rope;
		Rope hi_res;
		append(hi_res, lr);
		append(hi_res, rr);
		res.hi_rope = hi_res;
                res.hi_i_j_rope = le.hi_i_j_rope + re.hi_i_j_rope;

		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;

		string o1;
		string o2;
		string o3;
		append(o1, "acomb{", 6);
		//append(o1, res.firstStem);
		//append(o1, res.subword);
		//append(o1, res.lastStem);
		append(o1, le.rep);
		append(o2, ",", 1);
		append(o2, re.rep);
		append(o2, ",min(dr_energy(", 15);
		append(o2, dr_energy(le.lastStem, le.lastStem));
		append(o3, "),dl_energy(", 12);
		append(o3, dl_energy(re.firstStem, re.firstStem));
		append(o3, ")=", 2);
		append(o3, min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem)));
		append(o3, "}", 1);
		string o;
		append(o, o1);
		append(o, o2);
		append(o, o3);
                append(o, "[", 1);
                append(o, o6);
                append(o, "]", 1);
		res.rep = o;

		return res;
	}

	choice [mfeanswer_v4] h([mfeanswer_v4] i) {
		return list(minimum(i));
	}
}






algebra count auto count ;

algebra enum auto enum ;

algebra p_func implements Canonical_Algebra(alphabet = char, answer = pfanswer) {
	pfanswer sadd(Subsequence lb,pfanswer e) {
		pfanswer res = e;
		
		res.pf.q1 = scale(1) * e.pf.q1;
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer cadd(pfanswer le,pfanswer re) {
		pfanswer res = le;
		
		res.pf.q1 = le.pf.q1 * re.pf.q1;
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer cadd_Pr(pfanswer le,pfanswer re) {
		pfanswer res = le;
		
		res.pf.q1 = le.pf.q1 * sum_elems(re.pf);
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer cadd_Pr_Pr(pfanswer le,pfanswer re) {
		pfanswer res = le;
		
		res.pf = mk_tuple(le.firststem, le.pf.q1 * re.pf.q1);
		
		return res;
	}

	pfanswer cadd_Pr_Pr_Pr(pfanswer le,pfanswer re) {
		pfanswer res = le;
		
		res.pf = mk_tuple(le.firststem, le.pf.q1 * sum_elems(re.pf));
		
		return res;
	}

	pfanswer ambd(pfanswer le,Subsequence b,pfanswer re) {
		pfanswer res = le;
		
		res.pf.q1 = scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf);
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer ambd_Pr(pfanswer le,Subsequence b,pfanswer re) {
		pfanswer res = le;
		
		res.pf = mk_tuple(le.firststem, scale(1) * check_tuple(le.pf.q1, le.firststem, re.firststem, b, re.pf));
		
		return res;
	}

	pfanswer nil(Subsequence loc) {
		pfanswer res;
		
		res.pf.q1 = 1.0;
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		res.firststem.i = seq_size(loc);
		res.firststem.j = seq_size(loc);
		res.firststem.seq = loc.seq;
		
		return res;
	}

	pfanswer nil_Pr(Subsequence loc) {
		pfanswer res;
		
		res.pf.q1 = 1.0;
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		res.firststem.i = seq_size(loc);
		res.firststem.j = seq_size(loc);
		res.firststem.seq = loc.seq;
		
		return res;
	}

	pfanswer edl(Subsequence lb,pfanswer e) {
		pfanswer res = e;
		
		res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dl_energy(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer edr(pfanswer e,Subsequence rb) {
		pfanswer res = e;
		
		res.pf.q1 = scale(1) * e.pf.q1 * mk_pf(dr_energy(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer edlr(Subsequence lb,pfanswer e,Subsequence rb) {
		pfanswer res = e;
		
		res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(ext_mismatch_energy(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer drem(pfanswer e) {
		return e;
	}

	pfanswer is(pfanswer e) {
		pfanswer res = e;
		
		res.pf.q1 = e.pf.q1 * mk_pf(termau_energy(e.firststem, e.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer sr(Subsequence lb,pfanswer e,Subsequence rb) {
		pfanswer res = e;
		
		res.firststem.i = lb.i;
		res.firststem.j = rb.j;
		
		res.pf.q1 = scale(2) * e.pf.q1 * mk_pf(sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		pfanswer res;
		
		res.firststem.seq = llb.seq;
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(region.j - region.i + 4) * mk_pf(hl_energy(region) + sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer sp(Subsequence llb,Subsequence lb,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		res.pf.q1 = scale(4) * e.pf.q1 * mk_pf(sr_energy(res.firststem,res.firststem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer bl(Subsequence lregion,pfanswer e) {
		pfanswer res = e;
		
		res.firststem.i = lregion.i;
		
		Subsequence innerstem;
		innerstem.seq = lregion.seq;
		innerstem.i = lregion.i-1;
		innerstem.j = e.firststem.j+1;
		
		res.pf.q1 = scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(bl_energy(lregion,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer br(pfanswer e,Subsequence rregion) {
		pfanswer res = e;
		
		res.firststem.j = rregion.j;
		
		Subsequence innerstem;
		innerstem.seq = rregion.seq;
		innerstem.i = e.firststem.i-1;
		innerstem.j = rregion.j+1;
		
		res.pf.q1 = scale(rregion.j - rregion.i) * e.pf.q1 * mk_pf(br_energy(innerstem, rregion));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;

		return res;
	}

	pfanswer il(Subsequence lregion,pfanswer e,Subsequence rregion) {
		pfanswer res = e;
		
		res.firststem.i = lregion.i;
		res.firststem.j = rregion.j;
		
		res.pf.q1 = scale((lregion.j - lregion.i) + (rregion.j - rregion.i)) * e.pf.q1 * mk_pf(il_energy(lregion, rregion));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer ml(Subsequence llb,Subsequence lb,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(4) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}
	
	pfanswer mldr(Subsequence llb,Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(5) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dri_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladr(Subsequence llb,Subsequence lb,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q3) * mk_pf(min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(5) * amdangle * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldlr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(6) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladlr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		float amdangle;
		amdangle = e.pf.q1 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q2 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q3 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg( wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   e.pf.q4 * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)) + min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem, res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldladr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t rightdanglingBase = base_t(dr[dr.i]);
		base_t rightmostBaselastStem = base_t(e.firststem[dr.i-1]);
		double amdangle;
		amdangle = (e.pf.q1 * mk_pf(dli_energy(innerstem,innerstem)) + e.pf.q3) * mk_pf(min(dr_dangle_dg(wc_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem))) +
			   (e.pf.q2 + e.pf.q4) * mk_pf(min(dr_dangle_dg(wob_comp(rightmostBaselastStem), rightmostBaselastStem, rightdanglingBase), dri_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladldr(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem))) +
			   (e.pf.q3 + e.pf.q4 * mk_pf(dri_energy(innerstem,innerstem))) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(6) * amdangle * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mldl(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.pf.q1 = scale(5) * sum_elems(e.pf) * mk_pf(ml_energy() + ul_energy() + dli_energy(innerstem,innerstem) + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer mladl(Subsequence llb,Subsequence lb,Subsequence dl,pfanswer e,Subsequence rb,Subsequence rrb) {
		pfanswer res = e;
		
		res.firststem.i = llb.i;
		res.firststem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		base_t leftdanglingBase = base_t(dl[dl.i]);
		base_t leftmostBasefirstStem = base_t(e.firststem[dl.i+1]);
		float amdangle;
		amdangle = (e.pf.q1 + e.pf.q2) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem,  wc_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem))) +
			   (e.pf.q3 + e.pf.q4) * mk_pf(min(dl_dangle_dg(leftdanglingBase, leftmostBasefirstStem, wob_comp(leftmostBasefirstStem)), dli_energy(innerstem,innerstem)));
		
		res.pf.q1 = scale(5) * amdangle * mk_pf(ml_energy() + ul_energy() + sr_energy(res.firststem,res.firststem) + termau_energy(innerstem,innerstem));
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer addss(pfanswer e,Subsequence rregion) {
		pfanswer res = e;
		
		res.pf = mult_tup(scale(rregion.j - rregion.i) * mk_pf(ss_energy(rregion)), e.pf);

		return res;
	}

	pfanswer ssadd(Subsequence lregion,pfanswer e) {
		pfanswer res = e;
		
		Subsequence test;
		test.seq = lregion.seq;
		test.i = lregion.i;
		test.j = lregion.j+1;

		res.pf = mk_tuple(e.firststem, scale(lregion.j - lregion.i) * e.pf.q1 * mk_pf(ul_energy() + ss_energy(lregion)));
		
		return res;
	}

	pfanswer trafo(pfanswer e) {
		pfanswer res = e;
		
		res.pf.q1 = sum_elems(e.pf);
		res.pf.q2 = 0.0;
		res.pf.q3 = 0.0;
		res.pf.q4 = 0.0;
		
		return res;
	}

	pfanswer incl(pfanswer e) {
		pfanswer res = e;
		
		res.pf = mk_tuple(e.firststem, e.pf.q1 * mk_pf(ul_energy()));

		return res;
	}

	pfanswer combine(pfanswer le,pfanswer re) {
		pfanswer res = le;
		
		res.firststem = le.firststem;
		
		res.pf.q1 = (le.pf.q1 + le.pf.q2) * (re.pf.q1 + re.pf.q3);
		res.pf.q2 = (le.pf.q1 + le.pf.q2) * (re.pf.q2 + re.pf.q4);
		res.pf.q3 = (le.pf.q3 + le.pf.q4) * (re.pf.q3 + re.pf.q1);
		res.pf.q4 = (le.pf.q4 + le.pf.q3) * (re.pf.q4 + re.pf.q2);
		
		return res;
	}

	pfanswer acomb(pfanswer le,Subsequence b,pfanswer re) {
		pfanswer res = le;
		
		res.firststem = le.firststem;
		
		base_t baseLeftStem = base_t(le.firststem[b.i-1]);
		base_t baseRightStem = base_t(re.firststem[b.i+1]);
		base_t baseAmbigious = base_t(b[b.i]);
		double  wcDr = dr_dangle_dg(  wc_comp(baseLeftStem), baseLeftStem, baseAmbigious);
		double wobDr = dr_dangle_dg( wob_comp(baseLeftStem), baseLeftStem, baseAmbigious);
		double  wcDl = dl_dangle_dg(baseAmbigious, baseRightStem,  wc_comp(baseRightStem));
		double wobDl = dl_dangle_dg(baseAmbigious, baseRightStem, wob_comp(baseRightStem));
		
		res.pf.q1 = le.pf.q1 * (re.pf.q1 * mk_pf(min( wcDr, wcDl)) + re.pf.q3 * mk_pf(min( wcDr,wobDl))) + 
			    le.pf.q2 * (re.pf.q1 * mk_pf(min(wobDr, wcDl)) + re.pf.q3 * mk_pf(min(wobDr,wobDl)));
		res.pf.q2 = le.pf.q2 * (re.pf.q2 * mk_pf(min(wobDr, wcDl)) + re.pf.q4 * mk_pf(min(wobDr,wobDl))) + 
			    le.pf.q1 * (re.pf.q2 * mk_pf(min( wcDr, wcDl)) + re.pf.q4 * mk_pf(min( wcDr,wobDl)));
		res.pf.q3 = le.pf.q3 * (re.pf.q3 * mk_pf(min( wcDr,wobDl)) + re.pf.q1 * mk_pf(min( wcDr, wcDl))) +
			    le.pf.q4 * (re.pf.q3 * mk_pf(min(wobDr,wobDl)) + re.pf.q1 * mk_pf(min(wobDr, wcDl)));
		res.pf.q4 = le.pf.q4 * (re.pf.q4 * mk_pf(min(wobDr,wobDl)) + re.pf.q2 * mk_pf(min(wobDr, wcDl))) +
			    le.pf.q3 * (re.pf.q4 * mk_pf(min( wcDr,wobDl)) + re.pf.q2 * mk_pf(min( wcDr, wcDl)));
		
		res.pf.q1 = res.pf.q1 * scale(1);
		res.pf.q2 = res.pf.q2 * scale(1);
		res.pf.q3 = res.pf.q3 * scale(1);
		res.pf.q4 = res.pf.q4 * scale(1);
		
		return res;
	}

	choice [pfanswer] h([pfanswer] i) {
		return list(sum(i));
		//~ return i;
	}
}

algebra p_func_filter_me extends p_func {
  choice [pfanswer] h([pfanswer] l)
  {
    return l;
  }
}

//algebra p_func_sample extends p_func {
//  scoring choice [double] h([double] l)
//  {
//    return list(sample_value(l));
//  }
//}
//algebra p_func_id extends p_func {
//  choice [double] h([double] l)
//  {
//    return l;
//  }
//}

algebra p_func_id extends p_func {
  choice [pfanswer] h([pfanswer] l)
  {
    return l;
  }
}







algebra mfe implements Canonical_Algebra(alphabet = char, answer = mfeanswer) {
	mfeanswer sadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = e.energy;
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = e.firstStem.j;
		return res;
	}

	mfeanswer cadd(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer cadd_Pr_Pr_Pr(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer ambd_Pr(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		return res;
	}

	mfeanswer nil(Subsequence loc) {
		mfeanswer res;
		res.energy = 0;
		res.firstStem = loc;
		return res;
	}

	mfeanswer nil_Pr(Subsequence loc) {
		mfeanswer res;
		res.energy = 0;
		res.firstStem = loc;
		return res;
	}

	mfeanswer edl(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = e.energy + dl_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edr(mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + dr_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer edlr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + ext_mismatch_energy(e.firstStem, e.firstStem);
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer drem(mfeanswer e) {
		return e;
	}

	mfeanswer is(mfeanswer e) {
		mfeanswer res;
		res.energy = e.energy + termau_energy(e.firstStem, e.firstStem);
                //res.energy = e.energy;
		res.firstStem = e.firstStem;
		return res;
	}

	mfeanswer sr(Subsequence lb,mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = hl_energy(region) + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer sp(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		return res;
	}

	mfeanswer bl(Subsequence lregion,mfeanswer e) {
		mfeanswer res;
		res.firstStem.seq = lregion.seq;
		res.firstStem.i = lregion.i;
		res.firstStem.j = e.firstStem.j;
		
		Subsequence innerstem;
		innerstem.seq = lregion.seq;
		innerstem.i = lregion.i-1;
		innerstem.j = e.firstStem.j+1;
		
		res.energy = e.energy + bl_energy(lregion, innerstem);
		return res;
	}

	mfeanswer br(mfeanswer e,Subsequence rregion) {
		mfeanswer res;
		res.firstStem.seq = rregion.seq;
		res.firstStem.i = e.firstStem.i;
		res.firstStem.j = rregion.j;
		
		Subsequence innerstem;
		innerstem.seq = rregion.seq;
		innerstem.i = e.firstStem.i-1;
		innerstem.j = rregion.j+1;
		
		res.energy = e.energy + br_energy(innerstem, rregion);  
		return res;
	}

	mfeanswer il(Subsequence lregion,mfeanswer e,Subsequence rregion) {
		mfeanswer res;
		res.firstStem.seq = lregion.seq;
		res.firstStem.i = lregion.i;
		res.firstStem.j = rregion.j;
		
		res.energy = e.energy + il_energy(lregion, rregion);
		return res;
	}

	mfeanswer ml(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mldr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mladr(Subsequence llb,Subsequence lb,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem, e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mldladr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + min(dri_energy(innerstem,innerstem), dr_energy(e.lastStem,e.lastStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mladldr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer e,Subsequence rb,Subsequence rrb) {
		mfeanswer res;
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + min(dli_energy(innerstem,innerstem), dl_energy(e.firstStem, e.firstStem)) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		return res;
	}

	mfeanswer addss(mfeanswer e,Subsequence rb) {
		mfeanswer res;
		res.energy = e.energy + ss_energy(rb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.lastStem;
		return res;
	}

	mfeanswer ssadd(Subsequence lb,mfeanswer e) {
		mfeanswer res;
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer trafo(mfeanswer e) {
		return e;
	}

	mfeanswer incl(mfeanswer e) {
		mfeanswer res;
		res.energy = ul_energy() + e.energy;
		
		res.firstStem = e.firstStem;
		res.lastStem = e.firstStem;
		return res;
	}

	mfeanswer combine(mfeanswer le,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy;
		
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	mfeanswer acomb(mfeanswer le,Subsequence b,mfeanswer re) {
		mfeanswer res;
		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.firstStem = le.firstStem;
		res.lastStem = re.lastStem;
		return res;
	}

	choice [mfeanswer] h([mfeanswer] i) {
		return list(minimum(i));
	}
}

algebra mfe_filter_me extends mfe {
    choice [mfeanswer] h([mfeanswer] l) {
	return l;
    }
}

algebra mfe_id extends mfe {
    choice [mfeanswer] h([mfeanswer] l) {
	return l;
    }
}

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

	string hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((",2);
		append(res, '.', size(region));
		append(res, "))",2);
		return res;
	}

	string sp(Subsequence llb,Subsequence lb,string e,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((",2);
		append(res, e);
		append(res, "))",2);
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

	string ml(Subsequence llb,Subsequence lb,string e,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, e);
		append(res, "))", 2);
		return res;
	}

	string mldr(Subsequence llb,Subsequence lb,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
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

	string mldlr(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, '.');
		append(res, "))", 2);
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

	string mldl(Subsequence llb,Subsequence lb,Subsequence dl,string e,Subsequence rb,Subsequence rrb) {
		string res;
		append(res, "((", 2);
		append(res, '.');
		append(res, e);
		append(res, "))", 2);
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
		return i;
	}
}




// #########################################################################################
// ########################### hishape*V4 ######################################################

algebra hishapehV4 implements Canonical_Algebra(alphabet = char, answer = hishapeanswer_v4) {
	hishapeanswer_v4 sadd(Subsequence lb,hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;

		Rope empty_shape;
		if (e.hi_rope == empty_shape ) {
			res.hi_rope = "_";
		}

		string o6;
		append(o6, '.');
		append(o6, e.pp);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 cadd(hishapeanswer_v4 le,hishapeanswer_v4 re) {
		hishapeanswer_v4 res = le;    // for firstStem and other default values

		if (re.hi_rope == "_" ) {
			res.hi_rope = le.hi_rope;
		} else {
			res.hi_rope = le.hi_rope + re.hi_rope;
		}

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 cadd_Pr(hishapeanswer_v4 le,hishapeanswer_v4 re) {
                hishapeanswer_v4 res = le;

                res.hi_rope = le.hi_rope + re.hi_rope;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 cadd_Pr_Pr(hishapeanswer_v4 le,hishapeanswer_v4 re) {
		hishapeanswer_v4 res = le;


		if (re.hi_rope == "_" ) {
			res.hi_rope = le.hi_rope;
		} else {
			res.hi_rope = le.hi_rope + re.hi_rope;
		}

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 cadd_Pr_Pr_Pr(hishapeanswer_v4 le,hishapeanswer_v4 re) {
                hishapeanswer_v4 res = le;

                res.hi_rope = le.hi_rope + re.hi_rope;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 ambd(hishapeanswer_v4 le,Subsequence b,hishapeanswer_v4 re) {
                hishapeanswer_v4 res = le;

                res.hi_rope = le.hi_rope + re.hi_rope;

		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 ambd_Pr(hishapeanswer_v4 le,Subsequence b,hishapeanswer_v4 re) {
                hishapeanswer_v4 res = le;    

                res.hi_rope = le.hi_rope + re.hi_rope;

		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 nil(Subsequence loc) {
		hishapeanswer_v4 res;
		
		Rope r;                     // r is empty
		res.hi_rope = r;

		string o6;
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 nil_Pr(Subsequence loc) {
		hishapeanswer_v4 res;
		
		Rope r;
		res.hi_rope = r;

		string o6;
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 edl(Subsequence lb,hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;

		string o6;
		append(o6, '.');
		append(o6, e.pp);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 edr(hishapeanswer_v4 e,Subsequence rb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, e.pp);
		append(o6, '.');
		res.pp = o6;

		

		return res;
	}

	hishapeanswer_v4 edlr(Subsequence lb,hishapeanswer_v4 e,Subsequence rb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 drem(hishapeanswer_v4 e) {
                hishapeanswer_v4 res = e;


		return res;
	}

	hishapeanswer_v4 is(hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;
		

		
		return res;
	}

	hishapeanswer_v4 sr(Subsequence lb,hishapeanswer_v4 e,Subsequence rb) {
		hishapeanswer_v4 res = e;
		



		string o6;
		append(o6, '(');
		append(o6, e.pp);
		append(o6, ')');
		res.pp = o6;



		return res;
	}

	hishapeanswer_v4 hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res;
		
		Rope hi_res;
		int pos;
		pos = (lb.i+rb.j+1)/2;
		if ( pos*2 > lb.i+rb.j+1 ) pos = pos - 1;  
		append(hi_res, pos);
		if ( pos*2 != lb.i+rb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
		append(hi_res, ',');
		res.hi_rope = hi_res;






		string o6;
		append(o6, "((",2);
		append(o6, '.', size(region));
		append(o6, "))",2);
		res.pp = o6;



		return res;
	}

	hishapeanswer_v4 sp(Subsequence llb,Subsequence lb,hishapeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		

		



		string o6;
		append(o6, "((",2);
		append(o6, e.pp);
		append(o6, "))",2);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 bl(Subsequence lregion,hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;
		

		
		string o6;
		append(o6, '.', size(lregion));
		append(o6, e.pp);
		res.pp = o6;


		return res;
	}

	hishapeanswer_v4 br(hishapeanswer_v4 e,Subsequence rregion) {
		hishapeanswer_v4 res = e;
		

		
		string o6;
		append(o6, e.pp);
		append(o6, '.', size(rregion));
		res.pp = o6;


		return res;
	}

	hishapeanswer_v4 il(Subsequence lregion,hishapeanswer_v4 e,Subsequence rregion) {
		hishapeanswer_v4 res = e;
		

		
		string o6;
		append(o6, '.', size(lregion));
		append(o6, e.pp);
		append(o6, '.', size(rregion));
		res.pp = o6;


		return res;
	}

	hishapeanswer_v4 ml(Subsequence llb,Subsequence lb,hishapeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mldr(Subsequence llb,Subsequence lb,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mladr(Subsequence llb,Subsequence lb,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mldladr(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mladldr(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, '.');
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mldl(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,hishapeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, "((", 2);
		append(o6, '.');
		append(o6, e.pp);
		append(o6, "))", 2);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 addss(hishapeanswer_v4 e,Subsequence rb) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, e.pp);
		append(o6, '.', size(rb));
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 ssadd(Subsequence lb,hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;
		
		string o6;
		append(o6, '.', size(lb));
		append(o6, e.pp);
		res.pp = o6;
		

		return res;
	}

	hishapeanswer_v4 trafo(hishapeanswer_v4 e) {
                hishapeanswer_v4 res = e;


		return res;
	}

	hishapeanswer_v4 incl(hishapeanswer_v4 e) {
		hishapeanswer_v4 res = e;
		

		

		return res;
	}

	hishapeanswer_v4 combine(hishapeanswer_v4 le,hishapeanswer_v4 re) {
		hishapeanswer_v4 res = le;
		res.hi_rope = le.hi_rope + re.hi_rope;

		string o6;
		append(o6, le.pp);
		append(o6, re.pp);
		res.pp = o6;

		return res;
	}

	hishapeanswer_v4 acomb(hishapeanswer_v4 le,Subsequence b,hishapeanswer_v4 re) {
		hishapeanswer_v4 res = le;
		res.hi_rope = le.hi_rope + re.hi_rope;
		
		string o6;
		append(o6, le.pp);
		append(o6, '.');
		append(o6, re.pp);
		res.pp = o6;
		

		return res;
	}


	choice [hishapeanswer_v4] h([hishapeanswer_v4] i) {
		//return list(minimum(i));
		return unique(i);
	}
}




// #########################################################################################
// ########################### mfe*V4 ######################################################

algebra mfehV4 implements Canonical_Algebra(alphabet = char, answer = mfeanswer_v4) {
	mfeanswer_v4 sadd(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		Rope empty_shape;
                Rope underscore_shape;
                //append(underscore_shape, '_');
		//if (is_empty(e)) {
                //    append(res.hi_i_j_vector, underscore_shape);
		//}
		res.subword.i = lb.i;

		return res;
	}

	mfeanswer_v4 cadd(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		//if (is_empty(re)) {
                //    res.hi_i_j_vector = le.hi_i_j_vector;
		//} else {
                    res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;
		//}


		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 cadd_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
		//if (is_empty(re)) {
                //        res.hi_i_j_vector = le.hi_i_j_vector;
		//} else {
                        res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;
		//}


		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 cadd_Pr_Pr_Pr(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 ambd(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 ambd_Pr(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy + min(dr_energy(le.firstStem, le.firstStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;

		return res;
	}

	mfeanswer_v4 nil(Subsequence loc) {
		mfeanswer_v4 res;          // based on empty
		
		res.energy = 0;
		res.firstStem = loc;
		res.lastStem = loc;
		res.subword = loc;

		return res;
	}

	mfeanswer_v4 nil_Pr(Subsequence loc) {
		mfeanswer_v4 res;           // based on empty
		
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
		return e;
	}

	mfeanswer_v4 is(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
                if (!e._multiloop) {
                    Rope new_triplet;
		    append(new_triplet, e.triplet);
		    append(new_triplet, e.i1);
		    append(new_triplet, '%');
		    append(new_triplet, e.j1);
		    append(res.hi_i_j_vector, new_triplet);
                }

		res.energy = e.energy + termau_energy(e.firstStem, e.firstStem);  // for the AU pair additional energy
//res.empty_ = false;
		
		return res;
	}

	mfeanswer_v4 sr(Subsequence lb,mfeanswer_v4 e,Subsequence rb) {
		mfeanswer_v4 res = e;
		
                res.i1 = lb.i+1;
                res.j1 = rb.j;

		res.firstStem.seq = lb.seq;
		res.firstStem.i = lb.i;
		res.firstStem.j = rb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res;                // based on empty
		
                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, ',');
                res.triplet = hi_res;

                res.i1 = llb.i+1;
                res.j1 = rrb.j;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = hl_energy(region) + sr_energy(res.firstStem,res.firstStem);
		res.lastStem = res.firstStem;    // lastStem and subword have the same content of firstStem
		res.subword = res.firstStem;     // lastStem and subword have the same content of firstStem

		return res;
	}

	mfeanswer_v4 sp(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.i1 = llb.i+1;
                res.j1 = rrb.j;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
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
		
		res.energy = e.energy + bl_energy(lregion,innerstem);
		res.subword.i = lregion.i;
		
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
		
		return res;
	}

	mfeanswer_v4 il(Subsequence lregion,mfeanswer_v4 e,Subsequence rregion) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + il_energy(lregion, rregion);
		res.subword.i = lregion.i;
		res.subword.j = rregion.j;
		
		return res;
	}

	mfeanswer_v4 ml(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;
		
		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;
		
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

	mfeanswer_v4 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;
		
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

                res._multiloop = true;
		
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
		
                res._multiloop = true;

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

	mfeanswer_v4 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;
		
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

		return res;
	}

	mfeanswer_v4 ssadd(Subsequence lb,mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		
		res.energy = ul_energy() + e.energy + ss_energy(lb);
		res.subword.i = lb.i;
		
		return res;
	}

	mfeanswer_v4 trafo(mfeanswer_v4 e) {
		return e;
	}

	mfeanswer_v4 incl(mfeanswer_v4 e) {
		mfeanswer_v4 res = e;
		res.energy = ul_energy() + e.energy;  // TODO: what is ul energy?
		return res;
	}

	mfeanswer_v4 combine(mfeanswer_v4 le,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy;
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	mfeanswer_v4 acomb(mfeanswer_v4 le,Subsequence b,mfeanswer_v4 re) {
		mfeanswer_v4 res = le;
		
                res.hi_i_j_vector = le.hi_i_j_vector + re.hi_i_j_vector;

		res.energy = le.energy + re.energy + min(dr_energy(le.lastStem, le.lastStem), dl_energy(re.firstStem, re.firstStem));
		res.lastStem = re.lastStem;
		res.subword.j = re.subword.j;
		
		return res;
	}

	choice [mfeanswer_v4] h([mfeanswer_v4] i) {
		return list(minimum(i));
                //return i;
	}
}


algebra mfehplusV4 extends mfehV4 {
	mfeanswer_v4 ml(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

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

	mfeanswer_v4 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

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
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

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
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

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

	mfeanswer_v4 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res._multiloop = true;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

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
}

algebra mfemV4 extends mfehV4 {
	mfeanswer_v4 ml(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mldr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dri_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladr(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

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

	mfeanswer_v4 mldlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + ml_mismatch_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladlr(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;

                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

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
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

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
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

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

	mfeanswer_v4 mldl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		Subsequence innerstem;
		innerstem.seq = lb.seq;
		innerstem.i = lb.i;
		innerstem.j = rb.j;
		
		res.energy = ml_energy() + ul_energy() + e.energy + dli_energy(innerstem,innerstem) + sr_energy(res.firstStem,res.firstStem) + termau_energy(innerstem,innerstem);
		res.lastStem = res.firstStem;
		res.subword = res.firstStem;
		
		return res;
	}

	mfeanswer_v4 mladl(Subsequence llb,Subsequence lb,Subsequence dl,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                clear(res.hi_i_j_vector);
                append(res.hi_i_j_vector, "(");
                append(res.hi_i_j_vector, e.hi_i_j_vector);
                append(res.hi_i_j_vector, ")");

                Rope hi_res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "m,", 2);
                res.triplet = hi_res;

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
}

//                                 28b-1-55 26i-3-49 25b-9-41  27i-16-38 27-22-32  ==> 1-55 is not always the same as 1%55
//( ( 28b,26i,25b,27i,27, , (-350, 28b%3%49,26i%9%41,25b%16%38,27i%22%32,27%1%55,) ) ,    1
//12345678901234567890123456789012345678901234567890123456        
//(((((...(((....(((((.((((...)))).)))))))).....)))....)).
algebra mfebV4 extends mfemV4 {
	mfeanswer_v4 sp(Subsequence llb,Subsequence lb,mfeanswer_v4 e,Subsequence rb,Subsequence rrb) {
		mfeanswer_v4 res = e;
		
                res.i1 = llb.i+1;
                res.j1 = rrb.j;

//                Rope triplet_bi;
//		append(triplet_bi, e.triplet_bi);
//		append(triplet_bi, llb.i+1);
//		append(triplet_bi, '%');
//		append(triplet_bi, rrb.j);
//		append(res.hi_i_j_vector, triplet_bi);

		res.firstStem.seq = llb.seq;
		res.firstStem.i = llb.i;
		res.firstStem.j = rrb.j;
		
		res.energy = e.energy + sr_energy(res.firstStem,res.firstStem);
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
		
		res.energy = e.energy + bl_energy(lregion,innerstem);
		res.subword.i = lregion.i;
		
                Rope hi_res;
		int pos;
		pos = (res.subword.i+res.subword.j+1)/2;
		if ( pos*2 > res.subword.i+res.subword.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != res.subword.i+res.subword.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "b,", 2);
//                res.triplet_bi = hi_res;
                append(res.hi_i_j_vector, hi_res);  // added

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
		
                Rope hi_res;
		int pos;
		pos = (res.subword.i+res.subword.j+1)/2;
		if ( pos*2 > res.subword.i+res.subword.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != res.subword.i+res.subword.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "b,", 2);
//                res.triplet_bi = hi_res;
                append(res.hi_i_j_vector, hi_res);  // added

		return res;
	}

	mfeanswer_v4 il(Subsequence lregion,mfeanswer_v4 e,Subsequence rregion) {
		mfeanswer_v4 res = e;
		
		res.energy = e.energy + il_energy(lregion, rregion);
		res.subword.i = lregion.i;
		res.subword.j = rregion.j;

                Rope hi_res;
		int pos;
		pos = (res.subword.i+res.subword.j+1)/2;
		if ( pos*2 > res.subword.i+res.subword.j+1 ) pos = pos - 1;  
                append(hi_res, pos);
		if ( pos*2 != res.subword.i+res.subword.j+1 ) 
                {
                  append(hi_res, ".5", 2);
                }
                append(hi_res, "i,", 2);
//                res.triplet_bi = hi_res;
                append(res.hi_i_j_vector, hi_res);  // added
		
		return res;
	}
}

// #########################################################################################
// ########################### mfe*V4 end ##################################################




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

	Rope hl(Subsequence llb,Subsequence lb,Subsequence region,Subsequence rb,Subsequence rrb) {
		Rope res;
		//append(res, '[');
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
		append(res, ',');
		//append(res, ']');
		return res;
	}

	Rope sp(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
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

	Rope ml(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mldr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mladr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
		return e;
	}

	Rope mldlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
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

	Rope mldl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
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
	Rope ml(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
                Rope res;
                append(res, "(,", 2);
                append(res, e);
                append(res, "),", 2);
                return res;
	}

	Rope mldr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
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

	Rope mldlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
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

	Rope mldl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
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
	Rope ml(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
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

	Rope mldr(Subsequence llb,Subsequence lb,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
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

	Rope mldlr(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence dr,Subsequence rb,Subsequence rrb) {
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

	Rope mldl(Subsequence llb,Subsequence lb,Subsequence dl,Rope e,Subsequence rb,Subsequence rrb) {
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
	Rope sp(Subsequence llb,Subsequence lb,Rope e,Subsequence rb,Subsequence rrb) {
		Rope res;
		int pos;
		pos = (llb.i+rrb.j+1)/2;
		if ( pos*2 > llb.i+rrb.j+1 ) pos = pos - 1;  
		append(res, pos);
		if ( pos*2 != llb.i+rrb.j+1 ) append(res, ".5", 2);
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




grammar gra_microstate_rule1 uses Canonical_Algebra(axiom = struct) {

  struct    = cadd(dangle, struct) |     
              sadd(BASE, struct) |           
              nil(LOC) # h;

  dangle    = edl (BASE, initstem) suchthat mfe_thresh |                // rule 1  
              edr (initstem, BASE) suchthat mfe_thresh |                // rule 1  
              edlr(BASE, initstem, BASE) suchthat mfe_thresh |          // rule 1  
              drem(initstem) suchthat mfe_thresh # h;                   // rule 1  

  initstem = is(closed) # h;

  closed = stack | hairpin | leftB | rightB | iloop | multiloop # h;    // because of all the elements are ensured with stackpairing, in this step it doesn't need stackpairing any more

  stack = sr(BASE, closed, BASE) with basepairing # h;

  hairpin = hl(BASE, BASE, REGION with minsize(3), BASE, BASE) with stackpairing # h;

  leftB = sp(BASE, BASE, bl(REGION, closed), BASE, BASE) with stackpairing # h;

  rightB = sp(BASE, BASE, br(closed, REGION), BASE, BASE) with stackpairing # h;

  iloop = sp(BASE, BASE, il(REGION with maxsize(30), closed, REGION with maxsize(30)), BASE, BASE) with stackpairing # h;

  multiloop = ml   (BASE, BASE,                          ml_comps,                          BASE, BASE) with stackpairing |
              mldl (BASE, BASE, BASE,                    ml_comps,                          BASE, BASE) with stackpairing |
              mldr (BASE, BASE,                          ml_comps, BASE,                    BASE, BASE) with stackpairing |
              mldlr(BASE, BASE, BASE,                    ml_comps, BASE,                    BASE, BASE) with stackpairing # h;

  ml_comps  = cadd(incl(dangle), ml_comps1) |
              sadd(BASE, ml_comps) # h;          

  ml_comps1 = cadd(incl(dangle), ml_comps1) |
              sadd(BASE, ml_comps1) |
              incl(dangle) |
              addss(incl(dangle), REGION) # h;
}


grammar gra_microstate_rule1_rule2 uses Canonical_Algebra(axiom = struct) {
//
  struct    = cadd(dangle, struct) suchthat helix_repeat |     
              sadd(BASE, struct) suchthat helix_repeat |              // rule 2  
              nil(LOC) # h;    // difference 1: LOC vs EMPTY

// edl->edlr, edr->edlr, drem->edlr
  dangle    = edl (BASE, initstem) suchthat mfe_thresh |                // rule 1  
              edr (initstem, BASE) suchthat mfe_thresh |                // rule 1  
              edlr(BASE, initstem, BASE) suchthat mfe_thresh |          // rule 1  
              drem(initstem) suchthat mfe_thresh # h;                   // rule 1  

//##  dangle_ml    = edl (BASE, initstem) suchthat mfe_thresh_ml |
//##              edr (initstem, BASE) suchthat mfe_thresh_ml |
//##              edlr(BASE, initstem, BASE) suchthat mfe_thresh_ml |
//##              drem(initstem) suchthat mfe_thresh_ml # h;

//    dangle    = edlr(BASE, closed, LOC ) |
//                edlr(LOC,  closed, BASE) |
//                edlr(BASE, closed, BASE) |
//                edlr(LOC,  closed, LOC ) # h;


//	edanglel_ml = edl(BASE, initstem) # h;
//	edangler_ml = edr(initstem, BASE) # h;
//	edanglelr_ml = edlr(BASE, initstem, BASE) # h;
//	nodangle_ml = drem(initstem) # h;
//	initstem = is(closed) # h;

// example:   helene = f (REGION, closed, REGION) suchthat in_mfe_thresh ;
  initstem = is(closed) # h;             // difference 2: closed vs initstem

//  closed    = {stack                        |
//               hairpin                      |
//               leftB                        |
//               rightB                       |
//               iloop                        |
//               multiloop} with stackpairing # h;
  closed = stack | hairpin | leftB | rightB | iloop | multiloop # h;    // because of all the elements are ensured with stackpairing, in this step it doesn't need stackpairing any more

  stack = sr(BASE, closed, BASE) with basepairing # h;

  hairpin = hl(BASE, BASE, REGION with minsize(3), BASE, BASE) with stackpairing # h;

  leftB = sp(BASE, BASE, bl(REGION, closed), BASE, BASE) with stackpairing # h;

  rightB = sp(BASE, BASE, br(closed, REGION), BASE, BASE) with stackpairing # h;

  iloop = sp(BASE, BASE, il(REGION with maxsize(30), closed, REGION with maxsize(30)), BASE, BASE) with stackpairing # h;

  //multiloop = ml( BASE, BASE, ml_comps, BASE, BASE) # h;   // difference_4

  multiloop = ml   (BASE, BASE,                          ml_comps,                          BASE, BASE) with stackpairing |
              mldl (BASE, BASE, BASE,                    ml_comps,                          BASE, BASE) with stackpairing |
              mldr (BASE, BASE,                          ml_comps, BASE,                    BASE, BASE) with stackpairing |
              mldlr(BASE, BASE, BASE,                    ml_comps, BASE,                    BASE, BASE) with stackpairing # h;

//##  ml_comps  = sadd(BASE, ml_comps)        |
//##              cadd(incl(dangle_ml), ml_comps1) # h;
//##  ml_comps1 = sadd(BASE, ml_comps1)       |
//##              cadd(incl(dangle_ml), ml_comps1) |
//##              incl(dangle_ml)                  |
//##              addss(incl(dangle_ml), REGION) # h;  // difference 3  incl == ul

// following rule ensures a multiloop contains at least two dangles
  ml_comps  = cadd(incl(dangle), ml_comps1) suchthat helix_repeat |  // rule 2
              sadd(BASE, ml_comps) suchthat helix_repeat # h;          

  ml_comps1 = cadd(incl(dangle), ml_comps1) suchthat helix_repeat |  // rule 2
              sadd(BASE, ml_comps1) suchthat helix_repeat |
              incl(dangle) suchthat helix_repeat |
              addss(incl(dangle), REGION) suchthat helix_repeat # h;  // difference 3  incl == ul_energy. because all have ul_energy, it doesn't matter
// ssadd should be expensive O(n4)?? Why??? using sadd(sadd(sadd...)) 
}




grammar gra_microstate2 uses Canonical_Algebra(axiom = struct) {
  struct    = sadd(BASE, struct)   |
              cadd(dangle, struct) |
              nil(LOC)           # h;

  dangle    = edl (BASE, closed, LOC ) |
              edr (LOC,  closed, BASE) |
              edlr(BASE, closed, BASE) |
              drem(LOC,  closed, LOC ) # h;

  closed    = stack                        |
              hairpin                      |
              leftB                        |
              rightB                       |
              iloop                        |
              multiloop # h;  //TODO: comparing with canonicals_nonamb: there doesn't use stackpairing   ==> DELETE with stackpairing

  stack     = sr   (BASE,                                closed,                                  BASE) with stackpairing # h;
  hairpin   = hl   (BASE, BASE,                          REGION with minsize(3),            BASE, BASE) with stackpairing # h;
  leftB     = sp(BASE, BASE, bl(REGION, closed), BASE, BASE)                                            with stackpairing # h;
  rightB    = sp(BASE, BASE, br(closed, REGION), BASE, BASE)                                            with stackpairing # h;
  //leftB     = bl   (BASE, BASE, REGION,                  closed,                            BASE, BASE) with stackpairing # h;
  //rightB    = br   (BASE, BASE,                          closed,   REGION,                  BASE, BASE) with stackpairing # h;
  iloop     = il   (BASE, BASE, REGION with maxsize(30), closed,   REGION with maxsize(30), BASE, BASE) with stackpairing # h;  // !!!!WRONG il(3 args)

  multiloop = ml   (BASE, BASE,                          ml_comps,                          BASE, BASE) with stackpairing |
              mldl (BASE, BASE, BASE,                    ml_comps,                          BASE, BASE) with stackpairing |
              mldr (BASE, BASE,                          ml_comps, BASE,                    BASE, BASE) with stackpairing |
              mldlr(BASE, BASE, BASE,                    ml_comps, BASE,                    BASE, BASE) with stackpairing # h;

  ml_comps  = sadd(BASE, ml_comps)        |
              cadd(incl(dangle), ml_comps1) # h;

  ml_comps1 = sadd(BASE, ml_comps1)       |
              cadd(incl(dangle), ml_comps1) |
              incl(dangle)                  |
              addss(incl(dangle), REGION)   # h;
}


grammar canonicals_nonamb uses Canonical_Algebra(axiom = struct) {
	struct = left_dangle | trafo(noleft_dangle) | left_unpaired # h;

	left_unpaired = sadd(BASE, left_unpaired) | sadd(BASE, left_dangle) # h;

	left_dangle = ambd(edanglel, BASE, noleft_dangle) suchthat helix_repeat | 
                      cadd_Pr(edanglel, {noleft_dangle | nil_Pr(LOC)}) suchthat helix_repeat | 
                      cadd(edanglelr, {left_dangle | left_unpaired}) suchthat helix_repeat | 
                      nil(LOC) # h;                                                        // rule 2

	noleft_dangle = cadd_Pr_Pr(edangler, {left_dangle | left_unpaired}) suchthat helix_repeat | 
                        cadd_Pr_Pr_Pr(nodangle, {noleft_dangle | nil_Pr(LOC)}) suchthat helix_repeat | 
                        ambd_Pr(nodangle, BASE, noleft_dangle) suchthat helix_repeat # h;  // rule 2

// TODO: add choice function above_zero
	edanglel = edl(BASE, initstem) suchthat mfe_thresh # h;         // rule 1

	edangler = edr(initstem, BASE) suchthat mfe_thresh # h;         // rule 1

	edanglelr = edlr(BASE, initstem, BASE) suchthat mfe_thresh # h; // rule 1

	nodangle = drem(initstem) suchthat mfe_thresh # h;              // rule 1

//	edanglel_ml = edl(BASE, initstem) # h;

//	edangler_ml = edr(initstem, BASE) # h;

//	edanglelr_ml = edlr(BASE, initstem, BASE) # h;

//	nodangle_ml = drem(initstem) # h;

	initstem = is(closed) # h;

	closed = stack | hairpin | multiloop | leftB | rightB | iloop # h;

	multiloop = {mldl(BASE, BASE, BASE, ml_comps1, BASE, BASE) | 
		mladl(BASE, BASE, BASE, ml_comps2, BASE, BASE) | 
		mldr(BASE, BASE, ml_comps3, BASE, BASE, BASE) | 
		mladr(BASE, BASE, ml_comps2, BASE, BASE, BASE) | 
		mldlr(BASE, BASE, BASE, ml_comps4, BASE, BASE, BASE) | 
		mladlr(BASE, BASE, BASE, ml_comps2, BASE, BASE, BASE) | 
		mldladr(BASE, BASE, BASE, ml_comps1, BASE, BASE, BASE) | 
		mladldr(BASE, BASE, BASE, ml_comps3, BASE, BASE, BASE)  | 
		ml(BASE, BASE, ml_comps2, BASE, BASE)} with stackpairing # h;

	ml_comps1 = combine(block_dl, no_dl_no_ss_end) suchthat helix_repeat | combine(block_dlr, dl_or_ss_left_no_ss_end) suchthat helix_repeat | acomb(block_dl, BASE, no_dl_no_ss_end) suchthat helix_repeat # h;  // rule 2

	ml_comps2 = combine(incl(nodangle), no_dl_no_ss_end) suchthat helix_repeat | combine(incl(edangler), dl_or_ss_left_no_ss_end) suchthat helix_repeat | acomb(incl(nodangle), BASE, no_dl_no_ss_end) suchthat helix_repeat # h;  // rule 2

	ml_comps3 = combine(incl(edangler), dl_or_ss_left_ss_end) suchthat helix_repeat | combine(incl(nodangle), no_dl_ss_end) suchthat helix_repeat | acomb(incl(nodangle), BASE, no_dl_ss_end) suchthat helix_repeat # h;   // rule 2  

	ml_comps4 = combine(block_dl, no_dl_ss_end) suchthat helix_repeat | combine(block_dlr, dl_or_ss_left_ss_end) suchthat helix_repeat | acomb(block_dl, BASE, no_dl_ss_end) suchthat helix_repeat # h;   // rule 2


	block_dl = ssadd(REGION, edanglel) | incl(edanglel) # h;

	block_dlr = ssadd(REGION, edanglelr) | incl(edanglelr) # h;


	no_dl_no_ss_end = ml_comps2 | incl(nodangle) # h;

	dl_or_ss_left_no_ss_end = ml_comps1 | block_dl # h;

	no_dl_ss_end = ml_comps3 | incl(edangler) | addss(incl(edangler), REGION) # h;

	dl_or_ss_left_ss_end = ml_comps4 | block_dlr | addss(block_dlr, REGION) # h;




	stack = sr(BASE, closed, BASE) with basepairing # h;

	hairpin = hl(BASE, BASE, REGION with minsize(3), BASE, BASE) with stackpairing # h;

	leftB = sp(BASE, BASE, bl(REGION, closed), BASE, BASE) with stackpairing # h;

	rightB = sp(BASE, BASE, br(closed, REGION), BASE, BASE) with stackpairing # h;

	iloop = sp(BASE, BASE, il(REGION with maxsize(30), closed, REGION with maxsize(30)), BASE, BASE) with stackpairing # h;

}




//gapc hishapes0.gap -i hishapeh_mfe_pp_microstate -o test_hishapeh_mfe_pp_microstate.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapeh_mfe_pp_microstate = gra_microstate ( (hishapeh * mfeV4) * pretty ) ;
//gapc hishapes0.gap -i hishapehplus_mfe_pp_microstate -o test_hishapehplus_mfe_pp_macrostate.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapehplus_mfe_pp_microstate = gra_microstate ( (hishapehplus * mfeV4) * pretty ) ;
//gapc hishapes0.gap -i hishapem_mfe_pp_microstate -o test_hishapem_mfe_pp_macrostate.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapem_mfe_pp_microstate = gra_microstate ( (hishapem * mfeV4) * pretty ) ;
//gapc hishapes0.gap -i hishapeb_mfe_pp_microstate -o test_hishapeb_mfe_pp_macrostate.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapeb_mfe_pp_microstate = gra_microstate ( (hishapeb * mfeV4) * pretty ) ;


//gapc hishapes0.gap -i hishapeh_mfe_pp_macrostate -o test_hishapeh_mfe_pp_macrostate.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapeh_mfe_pp_macrostate = gra_microstate ( (hishapeh * mfe) * pretty ) ;
//## gapc hishapes0.gap -i hishapehV4_mfeV4_pp_microstate -o test_hishapehV4_mfeV4_pp_microstate.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapehV4_mfeV4_pp_microstate = gra_microstate ( (hishapehV4 * mfeV4) * pretty ) ;




//## gapc hishapes0.gap -i hishapehV4_mfeV4_pp -o hishapehV4_mfeV4_pp.cc -t --kbest --kbacktrace --no-coopt-class
//## make -f hishapehV4_mfeV4_pp.mf
//## ./hishapehV4_mfeV4_pp -f ../examples/collosoma_slrna.seq
instance hishapehV4_mfeV4_pp = canonicals_nonamb ( ((hishapehV4 * mfeV4) suchthat hishape_filter) * pretty); 


//## gapc hishapes0.gap -i hishapeh_extension_mfeV4_pp -o hishapeh_extension_mfeV4_pp.cc -t --kbest --kbacktrace --no-coopt-class
//## make -f hishapeh_extension_mfeV4_pp.mf
//## ./hishapeh_extension_mfeV4_pp -f ../examples/collosoma_slrna.seq
//D instance hishapeh_extension_mfeV4_pp = canonicals_nonamb ((hishapeh_extension * mfeV4) * pretty); 




// COMBINES FILTER METHOD 1 and the instance above
//## gapc hishapes0.gap -i hishapehV2_mfeV2_pp_filter_2 -o hishapehV2_mfeV2_pp_filter_2.cc -t --kbest
//## make -f hishapehV2_mfeV2_pp_filter_2.mf 
//instance hishapehV2_mfeV2_pp_filter_2 = canonicals_nonamb (hishapehV2 * (mfeV2 suchthat hishape_filter_5) );
instance hishapehV2_mfeV2_pp_filter_2 = canonicals_nonamb ( ((hishapeh * ( mfe % p_func )) suchthat pf_filter ) * pretty );


instance shapemfepfx = canonicals_nonamb ( ((hishapeh * ( mfe % p_func ))
            suchthat p_func_filter_all ) * pretty ) ;



// ##4## gapc hishapes0.gap -i hishapeh_mfe_pp -o hishapeh_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class 
// ##4## gapc hishapes0.gap -i hishapeh_mfe_pp -o rule1_h.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapeh_mfe_pp = gra_microstate ( (hishapeh * mfe) * pretty ) ;
// ##3## gapc hishapes0.gap -i hishapehplus_mfe_pp -o hishapehplus_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class
// ##3## gapc hishapes0.gap -i hishapehplus_mfe_pp -o rule1_hplus.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapehplus_mfe_pp = gra_microstate ( (hishapehplus * mfe) * pretty ) ;
// ##2## gapc hishapes0.gap -i hishapem_mfe_pp -o hishapem_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class
// ##2## gapc hishapes0.gap -i hishapem_mfe_pp -o rule1_m.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapem_mfe_pp = gra_microstate ( (hishapem * mfe) * pretty ) ;
// ##1## gapc hishapes.gap -i hishapeb_mfe_pp -o hishapeb_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class
// ##1## gapc hishapes.gap -i hishapeb_mfe_pp -o hishapeb_mfe_pp_for_lm.cc -t --kbacktrack --kbest --no-coopt-class  *
// ##1## gapc hishapes0.gap -i hishapeb_mfe_pp -o rule1_b.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapeb_mfe_pp = gra_microstate ( (hishapeb * mfe) * pretty ) ;




// ##4## gapc hishapes0.gap -i hishapeh_mfehV4_pp -o hishapeh_mfehV4_pp.cc -t --kbacktrack --kbest --no-coopt-class 
// ##4##for_test## gapc hishapes0.gap -i hishapeh_mfehV4_pp -o test_1.cc -t --kbacktrack --kbest --no-coopt-class 
// ##4##for_test## gapc hishapes0.gap -i hishapeh_mfehV4_pp -o test_5.cc -t --kbacktrack --kbest --no-coopt-class 
// ##4##for_test## gapc hishapes0.gap -i hishapeh_mfehV4_pp -o test_8.cc -t --kbacktrack --kbest --no-coopt-class 
// ##4##for_test## gapc hishapes.gap -i hishapeh_mfehV4_pp -o test_9.cc -t --kbacktrack --kbest --no-coopt-class 
// instance hishapeh_mfehV4_pp = gra_microstate ( (hishapeh * mfehV4) * pretty ) ;
// ##4##for_test## gapc hishapes0.gap -i hishapeh_mfehV4_pp_globalfilter -o test_3.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapeh_mfehV4_pp_globalfilter = gra_microstate ( (hishapeh * mfehV4) suchthat false_hishrep_filter * pretty ) ;


// new
// ##4## gapc hishapes.gap -i hishapeh_mfehV4_pp -o hishapeh_mfehV4_pp.cc -t --kbacktrack --kbest --no-coopt-class
// ##4## gapc hishapes.gap -i hishapeh_mfehV4_pp -o hishapeh_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class  
instance hishapeh_mfehV4_pp = gra_microstate ( (hishapeh * mfehV4) * pretty ) ;
// ##3## gapc hishapes.gap -i hishapehplus_mfehplusV4_pp -o hishapehplus_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapehplus_mfehplusV4_pp = gra_microstate ( (hishapehplus * mfehplusV4) * pretty ) ;
// ##2## gapc hishapes.gap -i hishapem_mfemV4_pp -o hishapem_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class 
instance hishapem_mfemV4_pp = gra_microstate ( (hishapem * mfemV4) * pretty ) ;
// ##1## gapc hishapes.gap -i hishapeb_mfebV4_pp -o hishapeb_mfebV4_pp.cc -t --kbacktrack --kbest --no-coopt-class
instance hishapeb_mfebV4_pp = gra_microstate ( (hishapeb * mfebV4) * pretty ) ;
// ##1## gapc hishapes.gap -i hishapeb_mfeV4_pp -o hishapeb_mfeV4_pp_for_lm.cc -t --kbacktrack --kbest --no-coopt-class  *
instance hishapeb_mfeV4_pp = gra_microstate ( (hishapeb * mfeV4) * pretty ) ;




// ##4## gapc hishapes.gap -i hishapeh_mfehV4_pfx -o hishapeh_mfe_pfx.cc -t --kbacktrack --kbest --no-coopt-class  
//instance hishapeh_mfehV4_pfx = canonicals_nonamb ( (hishapeh * mfehV4) * pretty ) ;
// ##3## gapc hishapes.gap -i hishapehplus_mfehplusV4_pfx -o hishapehplus_mfe_pfx.cc -t --kbacktrack --kbest --no-coopt-class
//instance hishapehplus_mfehplusV4_pfx = canonicals_nonamb ( (hishapehplus * mfehplusV4) * pretty ) ;
// ##2## gapc hishapes.gap -i hishapem_mfemV4_pfx -o hishapem_mfe_pfx.cc -t --kbacktrack --kbest --no-coopt-class 
//instance hishapem_mfemV4_pfx = canonicals_nonamb ( (hishapem * mfemV4) * pretty ) ;
// ##1## gapc hishapes.gap -i hishapeb_mfebV4_pfx -o hishapeb_mfe_pfx.cc -t --kbacktrack --kbest --no-coopt-class 
//instance hishapeb_mfebV4_pfx = canonicals_nonamb ( (hishapeb * mfebV4) * pretty ) ;

// new new
// ##4## gapc hishapes.gap -i hishapehV4_mfe_pp -o hishapehV4_mfe_pp.cc -t --kbacktrack --kbest --no-coopt-class  
instance hishapehV4_mfe_pp = canonicals_nonamb ((hishapehV4 * mfe) * pretty);


// gapc hix.gap -i hix_h_mfepfx -o hix_h_mfepfx.cc -t --kbacktrack --kbest ==> a lot faster (0m1.199s for k=100)
// ##4## gapc hishapes0.gap -i hishapeh_mfe_pfx -o hishapeh_mfe_pfx.cc -t --kbacktrack --kbest
instance hishapeh_mfe_pfx = canonicals_nonamb ( hishapeh *(mfe%p_func)*pretty ); 
// ##3## gapc hishapes0.gap -i hishapehplus_mfe_pfx -o hishapehplus_mfe_pfx.cc -t --kbacktrack --kbest
instance hishapehplus_mfe_pfx = canonicals_nonamb ( hishapehplus *(mfe%p_func)*pretty ) ;
// ##2## gapc hishapes0.gap -i hishapem_mfe_pfx -o hishapem_mfe_pfx.cc -t --kbacktrack --kbest
instance hishapem_mfe_pfx = canonicals_nonamb ( hishapem *(mfe%p_func)*pretty ) ;
// ##1## gapc hishapes0.gap -i hishapeb_mfe_pfx -o hishapeb_mfe_pfx.cc -t --kbacktrack --kbest
instance hishapeb_mfe_pfx = canonicals_nonamb ( hishapeb *(mfe%p_func)*pretty ) ;


instance p_func = canonicals_nonamb ( p_func );


// ##no kbest## gapc hishapes0.gap -i hishapem_mfe_pfx_nokbest -o hishapem_mfe_pfx_nokbest.cc -t --kbacktrack
instance hishapem_mfe_pfx_nokbest = canonicals_nonamb ( hishapem *(mfe%p_func)*pretty ) ;






// FILTER METHOD 2: doesn't work: using suchthat and no --kbest: gapc hishapes0.gap -i hishapeh_mfe_pfx_filter -o hishapeh_mfe_pfx_filter.cc -t --kbacktrack
//instance hishapeb_mfe_pfx_filter = canonicals_nonamb ( ((hishapeb *(mfe*p_func)) suchthat hishape_filter ) *pretty );    // consider hairpin, multi-, bulge- and internal-loops, additionally add a 'b' to end of bulge- and internal-loops
//instance hishapem_mfe_pfx_filter = canonicals_nonamb ( ((hishapem *(mfe%p_func)) suchthat hishape_filter ) *pretty );    // consider both hairpin- and multiloops, additionally add a 'm' to end of multiloops
instance hishapeh_mfe_pfx_filter = canonicals_nonamb ( (hishapeh *(mfe*p_func)  *pretty ) suchthat mfe_filter );    // only consider hairpin-loops
//1 instance hishapeh_mfe_pfx_filter = canonicals_nonamb ( (hishapeh suchthat hishape_filter)*(mfe%p_func)*pretty    );
//2 instance hishapeh_mfe_pfx_filter = canonicals_nonamb ( (hishapeh *(mfe%p_func)*pretty) suchthat hishape_filter  );

instance hishape_mfe_interleaved = canonicals_nonamb ( (hishapeh/mfe)*pretty);    // only consider hairpin-loops

// gapc -i hishapeh_mfe_pfx_filter -o hishapeh_mfe_pfx_filter.cc -t --kbacktrace --kbest hishapes0.gap  ?? TODO: check if --kbest and suchthat conflict
// make -f hishapemfepfx.mf
//
//instance hishapeh_mfe_pfx_filter = canonicals_nonamb ( ((hishapeh * ( mfe % p_func )) suchthat hishape_filter ) * pretty ) ;



// avoid use of hashtable: because of filterering #classes is
// relatively small -> constant factor of hashtable is
// significant in this usecase
//## gapc hishapes0.gap -i hishapeh_pfx_filter_using_hashtable -o hishapeh_pfx_filter_using_hashtable.cc -t --kbest
//## FILTER METHOD 1: make -f hishapeh_pfx_filter_using_hashtable.mf
instance hishapeh_pfx_filter_using_hashtable = canonicals_nonamb ((hishapeh * p_func_filter_me)
                                     suchthat hishape_filter_1);
