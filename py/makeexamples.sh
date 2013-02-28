## Automatic emaple generator 0.1 ( = m.tyka, procrastinating over his thesis ;) )

#frag/t_ff_gbsa.ppy      frag/t_head.ppy         frag/t_prot_minlangevin.ppy  frag/t_sys_aceglnnme.ppy  frag/t_sys_trpcage.ppy  frag/t_tail.ppy
#frag/t_ff_nb.ff         frag/t_prot_md.ppy      frag/t_prot_minsd.ppy        frag/t_sys_ala10.ppy      frag/t_sys_trpzip.ppy
#frag/t_ff_nb_dielec.ff  frag/t_prot_md_NVT.ppy  frag/t_prot_remd.ppy         frag/t_sys_protg.ppy      frag/t_sys_villin.ppy

cat frag/t_head.ppy frag/t_ff_gbsa.ppy       frag/t_sys_trpzip.ppy      frag/t_prot_minlangevin.ppy    frag/t_tail.ppy > langevin.py 
cat frag/t_head.ppy frag/t_ff_gbsa.ppy       frag/t_sys_villin.ppy      frag/t_prot_minlangevin.ppy    frag/t_tail.ppy > villlangevin.py 
cat frag/t_head.ppy frag/t_ff_gbsa.ppy       frag/t_sys_aceglnnme.ppy   frag/t_prot_remd.ppy           frag/t_tail.ppy > remd.py 
cat frag/t_head.ppy frag/t_ff_gbsa.ppy       frag/t_sys_aceglnnme.ppy   frag/t_prot_md_NVE.ppy         frag/t_tail.ppy > md2.py

cat frag/t_head.ppy frag/t_ff_nb_dielec.ppy  frag/t_sys_villin.ppy      frag/t_prot_md_NVE.ppy         frag/t_tail.ppy > md.py 

cat frag/t_head.ppy frag/t_ff_gbsa.ppy       frag/t_sys_allaa.ppy       frag/t_prot_epot.ppy           frag/t_tail.ppy > epot.py 
cat frag/t_head.ppy frag/t_ff_natcrest.ppy   frag/t_sys_trpcage.ppy     frag/t_prot_md_NVE.ppy         frag/t_tail.ppy > natc.py


cat frag/t_head.ppy frag/t_ff_nb_periodic.ppy frag/t_sys_tip3_15x15x15.ppy frag/t_prot_md_NVT.ppy   frag/t_tail.ppy > tip3box.py

 
