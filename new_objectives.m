% setting new objectives

% record indices of current objectives
ix_f = find(modelJoint.f);
ix_g = find(modelJoint.g);
ix_h = find(modelJoint.h);
ix_i = find(modelJoint.i);
ix_j = find(modelJoint.j);
ix_k = find(modelJoint.k);
ix_l = find(modelJoint.l);
ix_m = find(modelJoint.m);
ix_n = find(modelJoint.n);
ix_o = find(modelJoint.o);
ix_p = find(modelJoint.p);
ix_q = find(modelJoint.q);

% record indices for new objectives with rxn no.s
ix_new_f = 286;
ix_new_g = 1703;
ix_new_h = 3260;
ix_new_i = 5647;
ix_new_j = 6403;
ix_new_k = 7024;
ix_new_l = 7177;
ix_new_m = 8038;
ix_new_n = 12170;
ix_new_o = 13219;
ix_new_p = 15591;
ix_new_q = 22484;

% erase old indices to 0
modelJoint.f(ix_f) = 0;
modelJoint.g(ix_g) = 0;
modelJoint.h(ix_h) = 0;
modelJoint.i(ix_i) = 0;
modelJoint.j(ix_j) = 0;
modelJoint.k(ix_k) = 0;
modelJoint.l(ix_l) = 0;
modelJoint.m(ix_m) = 0;
modelJoint.n(ix_n) = 0;
modelJoint.o(ix_o) = 0;
modelJoint.p(ix_p) = 0;
modelJoint.q(ix_q) = 0;

% set new indices to 1
modelJoint.f(ix_new_f) = 1;
modelJoint.g(ix_new_g) = 1;
modelJoint.h(ix_new_h) = 1;
modelJoint.i(ix_new_i) = 1;
modelJoint.j(ix_new_j) = 1;
modelJoint.k(ix_new_k) = 1;
modelJoint.l(ix_new_l) = 1;
modelJoint.m(ix_new_m) = 1;
modelJoint.n(ix_new_n) = 1;
modelJoint.o(ix_new_o) = 1;
modelJoint.p(ix_new_p) = 1;
modelJoint.q(ix_new_q) = 1;