%|************************************************************************
%|SMARTI intersect
%|by Nils Reiners
%|Main source: Shirley, Marschner et al. 2010 - Fundamentals of computer graphics
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************


function [rst_struct hit]=intersect(p,P,S,rst_struct_init)

rst_struct=rst_struct_init;

count=0;
correct=0;
while count<3
    for n=1:length(S);
        
        P(3)=P(3)-correct;
        rst_struct(n).nr=n;
        
        v1=S(n).v1;
        v2=S(n).v2;
        p2=S(n).p;
        
        B=(p-p2);
        
        a=v1(1);           % just enumeration
        b=v1(2);
        c=v1(3);
        d=v2(1);
        e=v2(2);
        f=v2(3);
        g=-P(1);
        h=-P(2);
        i=-P(3);
        j=B(1);
        k=B(2);
        l=B(3);
        
        eimhf = e*i-h*f;    % calculate everything just once
        gfmdi = g*f-d*i;
        dhmeg = d*h-e*g;
        akmjb = a*k-j*b;
        jcmal = j*c-a*l;
        blmkc = b*l-k*c;
        
        M = a*eimhf + b*gfmdi + c*dhmeg;        % calculate determinant
        
        r = (j*eimhf+k*gfmdi+l*dhmeg)/M;     % reuse precalculations
        s = (i*akmjb+h*jcmal+g*blmkc)/M;
        t = -1*(f*akmjb+e*jcmal+d*blmkc)/M;
        
        
        rst_struct(n).r=r;
        rst_struct(n).s=s;
        rst_struct(n).t=t;
        
        if ~isfinite(rst_struct(n).r+rst_struct(n).s+rst_struct(n).t)==1
            rst_struct(n).r=0;
            rst_struct(n).s=0;
            rst_struct(n).t=0;
        end
        
        if rst_struct(n).t<=0.00001
            rst_struct(n).first_hit=0;
        elseif rst_struct(n).r>0&&rst_struct(n).s>0&&rst_struct(n).r+rst_struct(n).s<1;
            rst_struct(n).first_hit=1;
        else
            rst_struct(n).first_hit=0;
        end
        
    end
    
    %% Which surface?
    
    hit=find([rst_struct.first_hit]==1);
    
    if length(hit)>1
        
        N=find([rst_struct.first_hit]==1);
        O=rst_struct(N);
        Q=min([O.t]);
        hit=O(find([O.t]==Q)).nr;
        
    elseif isempty(hit)==1
        hit=0;
        return %Exits the function
    end
    if hit~=0
        break %Exits the while loop
    end
    correct=correct+0.0001;
    count=count+1;
end


