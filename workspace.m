vvv=repelem(G1,600);
rrr=repelem(G2,30);
 rrr=repmat(rrr,1,100);
 kkk=repmat(G3,1,2000);
 result=table(vvv',rrr',kkk',DIFC','VariableNames',{'valence' 'radius' 'koff' 'dif'});
 mod_result=result;
mod_result = sortrows(mod_result,"radius","ascend");
 %% 

dift=reshape(DIFC,30,2000 );
%x-valence ;y-radius ;different10 k off
figure()
mat_1=reshape(dift(1,:),20,100);
surf(G1,G2,mat_1)
xlabel("valence")
ylabel("np radius")
zlabel("prob cluster-uniform distribution")
legend()
title(['k off=0.1'])
figure()
mat_1=reshape(dift(15,:),20,100);
surf(G1,G2,mat_1)
xlabel("valence")
ylabel("np radius")
zlabel("prob cluster-uniform distribution")
legend()
title(['k off=2.8072'])
figure()
mat_1=reshape(dift(23,:),20,100);
surf(G1,G2,mat_1)
xlabel("valence")
ylabel("np radius")
zlabel("prob cluster-uniform distribution")
legend()
title(['k off=18.8739 '])
%% 
figure()
dift2=reshape(DIFC,600,100 )
mat_1=reshape(dift2(:,1),30,20)
surf(G2,G3,mat_1)
ylabel("k off")
xlabel("np radius")

zlabel("prob cluster-uniform distribution")
legend()
title(['valence=1'])
figure()
dift2=reshape(DIFC,600,100 )
mat_1=reshape(dift2(:,2),30,20)
surf(G2,G3,mat_1)
ylabel("k off")
xlabel("np radius")
zlabel("prob cluster-uniform distribution")
legend()
title(['valence=2'])
figure()
dift2=reshape(DIFC,600,100 )
mat_1=reshape(dift2(:,20),30,20)
surf(G2,G3,mat_1)
ylabel("k off")
xlabel("np radius")
zlabel("prob cluster-uniform distribution")
legend()
title(['valence=20'])
%% 
%x label valence y lable is k off 
figure()
dift3=reshape(table2array(mod_result(:,4)),3000,20)
mat_1=reshape(dift3(:,1),30,100)
surf(G1,G3,mat_1)
ylabel("k off")
xlabel("valence")
zlabel("prob cluster-uniform distribution")
legend()
title(['np_radius=50'])

figure()
dift3=reshape(table2array(mod_result(:,4)),3000,20)
mat_1=reshape(dift3(:,4),30,100)
surf(G1,G3,mat_1)
ylabel("k off")
xlabel("valence")
zlabel("prob cluster-uniform distribution")
legend()
title(['np radius=121.0526'])
figure()
dift3=reshape(table2array(mod_result(:,4)),3000,20)
mat_1=reshape(dift3(:,20),30,100)
surf(G1,G3,mat_1)
ylabel("k off")
xlabel("valence")
zlabel("prob cluster-uniform distribution")
legend()
title(['np radius=500'])
