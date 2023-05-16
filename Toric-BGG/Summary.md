# Differential modules and toric BGG.

Here is a sketch of our goals for the week. All the papers I'll reference are uploaded to our Github page.

# Minimal free resolutions of differential modules

The first goal will be to implement minimal free resolutions of differential modules, as introduced in the paper "Minimal free resolutions of differential modules" by Brown-Erman, which I'll abbreviate as BE1. We'll mainly focus on implementing

(a) Construction 2.7 in BE1, which is an algorithm for building a (possibly non-minimal) free resolution of a differential module.

(b) Construction 2.8 in BE1, which gives an alternative (and simpler) algorithm than (a).

(c) the ``degeneration" procedure in the proof of Theorem 3.2 in BE1.

(d) the proof of Proposition 4.1 in BE1, which is an algorithm for "minimizing" a non-minimal free differential module.

In fact, we may only end up implementing (b) and (d). For background on differential modules, I suggest skimming BE1, paying special attention to Construction 2.7 and the proof of Proposition 4.1.

# Toric BGG functors

Our second goal is to implement the toric BGG functors, as discussed in "Tate resolutions on toric varieties" by Brown-Erman (BE2). 

Specifically, we want to implement the functors L and R in Section 2.2 of BE2. Reading Section 2 of that paper up through Theorem 2.7 should give sufficient background. 

This has already been done in the standard graded case: see the BGG package, which is based on Eisenbud-Floystad-Schreyer's paper "Sheaf cohomology and free resolutions over exterior algebras". The idea is to extend many aspects of that package to the multigraded setting.

# Applications 

(1) Implement a method for computing "strongly linear strands" of multigraded free resolutions, as in the paper "Linear strands of multigraded free resolutions" by Brown-Erman (BE3). This should be really easy, once the code for the BGG functors is in place: we will just use the formula for the strongly linear strand in Theorem 1.3(2) of BE3. For background on this topic, reading the intro to BE3 should suffice.

(2) We will also implement Tate resolutions and sheaf cohomology computations over weighted projective stacks. This will require combining our minimal free resolutions of DM's code with the BGG code. The reference for this is Section 3 of BE2. But if you've never seen this material before, this reading may be heavy going, and so it may be best to just discuss this topic in person.
