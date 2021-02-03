                i=int(xt/dx+0.5)
                j=int(yt/dy+0.5)
                k=int(z2(ns,m)/dz+0.5)-gclr*kcnt

            exp1=exp1 + ex(i,j,k)
            eyp=eyp + ey(i,j,k)
            ezp =ezp + ez(i,j,k)
            delbxp = delbxp+delbx(i,j,k)     
            delbyp = delbyp+delby(i,j,k)     
            dpdzp = dpdzp+dpdz(i,j,k)     
            dadzp = dadzp+dadz(i,j,k)     
            aparp = aparp+apar(i,j,k)
            !twk: note the sign
            rhoreyp = rhoreyp - rhog*(e1rp*sn(l) + e2rp*cn(l))*ey(i,j,k)
            rhozeyp = rhozeyp - rhog*(e1zp*sn(l) + e2zp*cn(l))*ey(i,j,k)
            rhoztexp = rhoztexp - rhog*e2ztp*cn(l)*ex(i,j,k)
            rhozteyp = rhozteyp - rhog*e2ztp*cn(l)*ey(i,j,k)
            rhoztezp = rhoztezp - rhog*e2ztp*cn(l)*ez(i,j,k)
