function zcra=z(t,r,b_FR,p_FR,q_main)
zcra=((b_FR+2).*(1+(p_FR-1).*(r.^b_FR)).*q_main.*t./(b_FR+2.*p_FR));




