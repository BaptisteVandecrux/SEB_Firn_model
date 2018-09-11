function c = update_column_properties(c,psnowc, psnic, pslwc)

c.cdel=psnowc+pslwc+psnic;

c.cdelsum = 0;
aux = cumsum(c.cdel);

for jk = 1:c.jpgrnd
    c.cdelsum = c.cdelsum + c.cdel(jk);
    c.cmid(jk) = c.cdelsum - ( c.cdel(jk) / 2 );
    c.rcdel(jk) = 1/c.cdel(jk);
end
end