
function [EXP_numpChips,EXP_solvent,EXP_run,EXP_date,EXP_H_channel,EXP_W_channel,EXP_flowrate,EXP_desiredframerate,EXP_collection] = B2KExpParam(vid)
expression_numpChips = '(\d+)p_(\w+)_R(\d+)';
matchExp_numpChips = regexp(vid,expression_numpChips,'match');
EXP_numpChips = convertCharsToStrings(extractBefore(char(matchExp_numpChips),"p_"));

expression_solvent = 'p_(\w+)_R(\d+)';
matchExp_solvent = regexp(vid,expression_solvent,'match');
EXP_solvent = convertCharsToStrings(extractBetween(char(matchExp_solvent),"p_","_R"));

expression_run = '(\w+)_R(\d+)';
matchExp_run = regexp(vid,expression_run,'match');
EXP_run = str2double(extractAfter(char(matchExp_run),"_R"));

expression_date = 'R(\d+)_(\d+)_H(\d+).(\d+)_W(\d+).(\d+)';
matchExp_date = regexp(vid,expression_date,'match');
EXP_date = convertCharsToStrings(char(extractBetween(matchExp_date,"_","_")));

expression_channel = 'H(\d+).(\d+)_W(\d+).(\d+)';
matchExp_channel = regexp(vid,expression_channel,'match');
EXP_H_channel = str2double(extractBetween(char(matchExp_channel),"H","_W"));
EXP_W_channel = str2double(extractAfter(char(matchExp_channel),"_W"));

expression_flowrate = 'W(\d+).(\d+)_(\d+)mL-min';
matchExp_flowrate = regexp(vid,expression_flowrate,'match');
EXP_flowrate = str2double(extractBetween(char(matchExp_flowrate),"_","mL-min"));

expression_framerate = 'mL-min_(\d+)FPS';
matchExp_framerate = regexp(vid,expression_framerate,'match');
EXP_desiredframerate = str2double(extractBetween(char(matchExp_framerate),"mL-min_","FPS"));

expression_collection = 'FPS_(\d+)';
matchExp_collection = regexp(vid,expression_collection,'match');
EXP_collection = str2double(extractAfter(char(matchExp_collection),"FPS_"));
end
