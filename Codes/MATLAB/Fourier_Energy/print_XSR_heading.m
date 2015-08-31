%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Print Program Header
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function print_XSR_heading()

time = datevec(now);
time = get_time(time);

disp('***************************************************************************')
disp('***************************************************************************')
disp('*                                                                         *')
disp('*               ____________ _____    __   __ ___________                 *')
disp('*               | ___ \  _  \_   _|   \ \ / //  ___| ___ \                *')
disp('*               | |_/ / | | | | |______\ V / \ `--.| |_/ /                *')
disp('*               |  __/| | | | | |______/   \  `--. \    /                 *')
disp('*               | |   | |/ /  | |     / /^\ \/\__/ / |\ \                 *')
disp('*               \_|   |___/   \_/     \/   \/\____/\_| \_|                *')
disp('*                                                                         *')
disp('*                                                                         *')
disp('***************************************************************************')
disp('***************************************************************************')
disp('*                                                                         *')
disp('*                                                                         *')
disp('*  This code was created for the express use of author Michael Hackemack  *')
disp('*  and any and all reproduction, alteration or utilization of this source *')
disp('*  with the express use of redistribution in a commercial code without    *')
disp('*  prior approval of the author is in violation of human decency.         *')
disp('*                                                                         *')
disp('*  (c) 2015 - Michael Hackemack                                           *')
disp('*                                                                         *')
disp('***************************************************************************')
disp('                                                                           ')
disp(['  Current Run Date:  ',date])
disp(['  Current Run Time:  ',time])
disp('                                                                           ')
disp('---------------------------------------------------------------------------')
disp('======================== Begin Main Code Execution ========================')
disp('---------------------------------------------------------------------------')
disp('                                                                           ')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%	Function List                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_time(time)

hour = num2str(time(4));
minute = num2str(time(5));
second = num2str(round(time(6)));

if length(hour) == 1
    hour = ['0',hour];
end
if length(minute) == 1
    minute = ['0',minute];
end
if length(second) == 1
    second = ['0',second];
end

out = [hour,' : ',minute,' : ',second];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%