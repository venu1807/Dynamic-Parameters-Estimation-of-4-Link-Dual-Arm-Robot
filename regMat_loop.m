% This function evaluates the regressor matrix of a free floating robot
% according to the formulation in the paper. Please note that the second
% section has to be updated as per the robot's kinematic structure
% q, th_d, dq, dth_d - position and velocity of base and joints
% lk_no - index of link which is attached to the base

function [Pl] = regMat_loop(q,th_d,dq,dth_d,lk_no)
% Start of the actual program - ok
[n, ~, ~, ~, ~, ~, ~, ~, ~, ~, alt]=inputs();

% Update the position of sensor in base CM frame (sa0) - ok
% sa0 = [-0.07; -0.1; 0];
sa0 = [-0.2; -0.3; 0];
%Positions of the joints fixed to the base in sensor frame
% (in real system this will be an input)
jt_pos = [alt(1)-sa0(1), -sa0(1)-alt(1)
          -sa0(2)      , -sa0(2) 
          0            , 0             ];

% Robot's Kinematic Data (in inertial frame) - ok

q=load('statevar.dat');
t=load('timevar.dat');
m=load('mtvar.dat')/10;

k=length(t);
U=zeros(2*k,3*n);
V=zeros(2*k,1);

U1=zeros(2*k,(3*n)-1);
D=zeros(2*k,1);
count=1;

% q=[-1.624088e-09,	 3.776530e-09	, 0.000000e+00	, 2.918448e-08,	 0.000000e+00,	 0.000000e+00	, 3.802995e-01]	;
% dq=[-4.872189e-08,	 1.132943e-07,	 0.000000e+00,	 8.755216e-07,	 0.000000e+00,	 0.000000e+00,	 -1.426359e-05]	; 
 for t1=1:length(t)
    
    lk_no=[1];

    r0c = q(t1,1:3)'; %CM position of the base
    th0 = q(t1,4); %Angular position of the base
    th_d=q(t1,7);
    thj = th_d; %Joint position

    v0c = q(t1,8:10)'; %CM linear velocity of the base
    om0 = q(t1,11); %Angular velocity of the base
    dth_d=q(t1,14);
    omj = dth_d; %Joint velocity

    nc = length(lk_no); %number of chains attached to the base link

% Computation of sensor data - ok
    R=zeros(3,3,n);
    R(:,:,1) = [cos(th0) -sin(th0) 0
            sin(th0)  cos(th0) 0
            0         0        1];

% %Fabrication of the simulated sensor data
% (in real system this is the input instead of base CM pose and velocity)
    rs = r0c + R(:,:,1)*sa0;
    vs = v0c + cross([0;0;om0], R(:,:,1)*sa0); 
%

% Computation of the rotation matrices - ok
    for j=1:n-1
        R(:,:,j+1) = [cos(thj(j)) -sin(thj(j)) 0
                      sin(thj(j))  cos(thj(j)) 0
                      0              0         1];
    end

    RiI = zeros(3,3,n);
    RiI(:,:,1) = R(:,:,1);
    lkm = [lk_no, n];

    for j=1:nc
         for i=lkm(j):(lkm(j+1)-1)
           if i == lkm(j)
               mul = RiI(:,:,1);
           else
               mul = RiI(:,:,i);
           end
          RiI(:,:,i+1) = mul*R(:,:,i+1);
         end
    end

% Computation of joint vector(s) connected to the base in the inertial frame - ok
    rsiJ = zeros(3,nc); %Vector from sensor to joints fixed to the base

    for i=1:nc
             rsiJ(:,i) = RiI(:,:,1)*jt_pos(:,i);
    end

% Computation of joint position in the inertia frame - ok
    ri = zeros(3,n); % joint vectors from the origin of the inertial frame
    ri(:,1) = rs;

    count1 = 1;
    count2 = 2;
    for j=1:nc
       for i=lkm(j):(lkm(j+1)-1)
           if i == lkm(j)
                 add = rsiJ(:,count1);   
                 prev = ri(:,1);
                 count1 = count1+1;
           else
                 add = RiI(:,:,i)*[alt(count2);0;0];
                 prev = ri(:,i);
                 count2 = count2 + 1;
           end
          ri(:,i+1) = prev + add;
       end
        count2 = count2  + 1;
    end

% Computation of regressor matrix's elements
% vk and alpha computation - ok

    vk = zeros(3,n); 
    afa = zeros(3,n);
    vk(:,1) = vs;
    afa(:,1) = cross(ri(:,1),vk(:,1));
    exc = 0;
    for j=1:nc
         for i=lkm(j):(lkm(j+1)-1)
             vk(:,i+1) = vk(:,1) + cross([0;0;om0], (ri(:,i+1) - rs));
             if lkm(j) <= i-1
                  for h = lkm(j):i-1
                     rot = dth_d(h);
                     exc = exc + cross([0;0;rot], ri(:,i+1) - ri(:,h+1)-rs);
                  end
             end
             vk(:,i+1) = vk(:,i+1) + exc;
             afa(:,i+1) = cross(ri(:,i+1), vk(:,i+1));
             exc = 0;
         end
    end

% omk and beta computation

    omk = zeros(1,n);
    bomk = zeros(3,3,n);
    bta = zeros(3,n); 
    bbta = zeros(3,3,n);
    omk(1) = om0;
    bomk(:,:,1) = ([0,-omk(1),0;omk(1),0,0;0,0,1])*RiI(:,:,1);
% bta(:,1) = -vk(:,1) + cross(ri(:,1),[0;0;omk(1)]);
% bbta(:,:,1) = box(bta(:,1))*RiI(:,:,1);

    for j=1:nc
         for i=lkm(j):(lkm(j+1)-1)
              if i == lkm(j)
                  prev = omk(1);
              else
                  prev = omk(i);
              end
              omk(i+1) = prev + omj(i);
              bomk(:,:,i+1) = ([0,-omk(i+1),0;omk(i+1),0,0;0,0,1])*RiI(:,:,i+1);
%         bta(:,i+1) = -vk(:,i+1) + cross(ri(:,i+1),[0;0;omk(:,i+1)]);
%         bbta(:,:,i+1) = box(bta(:,i+1))*RiI(:,:,i+1);
         end
    end


% Deleting the zero components

    vk(3,:)=[]; 
    afa([1,2],:)=[];
    bomk(:,3,:) = [];
    bomk(3,:,:) = [];
% bbta([1,2],:,:) = []; 
% bbta(:,3,:) = [];


 % Matrix assembly for all the n links

    G1=zeros(2,2);
    GI=zeros(2,2);


    z21 = zeros(2,1);
    
    G1=[bomk(1,2,1) bomk(1,1,1);bomk(2,2,1) bomk(2,1,1)];
% prow = [vk(:,1)  bomk(:,:,1)];
    prow = [vk(:,1)  G1];
    lrow = omk(1);
    for i=2:n
         GI=[bomk(1,2,i) bomk(1,1,i);bomk(2,2,i) bomk(2,1,i)];
         papd = [vk(:,i)  GI];
    %  papd = [vk(:,i) bomk(:,:,i) ];
         prow = [prow, papd];
%     lrow = [lrow, omk(i)];
    end
     prow;
    
%      t2=1:3:601;
     U(count,:)=prow(1,:);
     U(count+1,:)=prow(2,:);
     
     V(count,1)  =m(t1,1);
     V(count+1,:)=m(t1,2);
     
     count=count+2;
 end
U;
V;
param=[500;500*0.3;500*0.2;10;10*0;10*0.5];
Pl=U*param;

for l=1:3*n
    if l<4
        U1(:,l)=U(:,l);
    elseif (l>4 && l<=(3*n))
        U1(:,l-1)=U(:,l);
    elseif l==4
        D=U(:,l);
    end
end

 U1;
 D;
    
    
    
% C=[U(:,1) U(:,2) U(:,3) U(:,5) U(:,6)];
% D=U(:,4);
% JJ=pinv(U1)*(V-D)

% for l=1:3*n
%      U1(:,l)=[U(:,l)]

 
% Verifiction of conservation of momentum
% ltum = prow*[500 0 0 10 5 0 10 5 0 10 5 0 10 5 0 10 5 0 10 5 0 ].';
% atum = low*[10.05 1.05 1.05].'
% prow
end
