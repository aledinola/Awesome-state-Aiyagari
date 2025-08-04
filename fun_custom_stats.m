function custom_stats = fun_custom_stats(V,Policy,StationaryDist,Params,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_grid,pi_z,heteroagentoptions,vfoptions,simoptions)
% CustomStats=Aiyagari1994_CustomModelStats(V,Policy,StationaryDist,Parameters,FnsToEvaluate,n_d,n_a,n_z,d_grid,a_grid,z_gridvals,pi_z,heteroagentoptions,vfoptions,simoptions)
% % CustomStats: output a structure with custom model stats by fields with the names you want
% % The inputs to CustomModelStats() are fixed and not something the user can
% % change. They differ based on InfHorz of FHorz, and differ if using PType.
% % Note: Most of them are familiar, only thing that may confuse is
% % 'z_gridvals' which is a joint-grid rather than whatever the user set up.
% CustomStats=struct();
% % Note: CustomStats must be scalar-valued

custom_stats = struct();

% corr_h_z
% cv_h

% --- TARGET 1 
AllStats=EvalFnOnAgentDist_AllStats_Case1(StationaryDist,Policy,FnsToEvaluate,...
    Params,[],n_d,n_a,n_z,d_grid,a_grid,z_grid,simoptions);

custom_stats.wealth_top1 = 1-AllStats.K.LorenzCurve(99);

end %end function