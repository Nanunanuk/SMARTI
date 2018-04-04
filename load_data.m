%|************************************************************************
%|SMARTI load_data
%|by Nils Reiners
%|Source: PV-lighthouse
%|This source code is free: you can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation., version 3 of the License. Any redistribution of		 |% %|modified code must also be free under the GNU General Public License.																																							 |%
%|This source code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for |% %|more details. You should have received a copy of the GNU General Public License along with this source code. If not, see <http://www.gnu.org/licenses/>.																		 |%
%|************************************************************************


function load_data(Mat)

%% Loading radiation data and material parameters
curr_path=mfilename('fullpath');
[a b c]=xlsread([curr_path(1:length(curr_path)-9) 'Data\All_data_PV-Lighthouse.xlsx']);


for n = 1:size(Mat,1)
    row1=find(strcmp(Mat{n,1},b(1,:)));
    row2=find(strcmp(Mat{n,2},b(2,:)));
    row3=find(strcmp(Mat{n,3},b(3,:)));
    
    x=row1(ismember(row1,row2));
    data_col=x(ismember(x,row3));
    
    wl_n_k=[a(:,1) a(:,data_col) a(:,data_col+1)];
    assignin('base', [Mat{n}], wl_n_k);
end