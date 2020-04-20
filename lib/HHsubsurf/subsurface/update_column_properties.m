function c = update_column_properties(c,psnowc, psnic, pslwc)

c.cdel=psnowc+pslwc+psnic;
aux = cumsum(c.cdel);

c.cdelsum = 0;

for jk = 1:c.jpgrnd
    c.cdelsum = c.cdelsum + c.cdel(jk);
    c.cmid(jk) = c.cdelsum - ( c.cdel(jk) / 2 );
    c.rcdel(jk) = 1/c.cdel(jk);
end
end