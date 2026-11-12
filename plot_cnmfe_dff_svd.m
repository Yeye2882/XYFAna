function [dff,F0,fig] = plot_cnmfe_dff_svd(neuron,Fs,varargin)
% 分块 SVD 自动版本：只需 plot_cnmfe_dff_svd(neuron,neuron.Fs)。
% 默认全时段背景中位数，与本例 df_prctile=50, df_window=[] 一致。
% 需要 P.mat_data 对应 MAT 文件中的 patch_pos；不需要读取原始电影。
% 不修改 neuron；结果写入输出变量 dff/F0。
if nargin<2 || isempty(Fs), Fs=neuron.Fs; end
%PLOT_CN MFE_DFF Normalize CNMF-E activity and plot a Delta F/F0 heatmap.
% 保存为 plot_cnmfe_dff_svd.m；先运行 cnmfe_setup，再 load 原始 MAT。
%
% 用法 (Fs 为活动矩阵的实际采样率，Hz):
%   [dff,F0,fig] = plot_cnmfe_dff_svd(neuron,10);
%   % 自动使用已有 neuron.C_df，或 neuron.Df。
%
%   [dff,F0,fig] = plot_cnmfe_dff_svd(neuron,10,'Ybg',Ybg);
%   % Ybg: CNMF-E 估计的完整背景电影，像素数 x 帧数；
%   % 也接受 高 x 宽 x 帧数，与 neuron.A 使用相同像素顺序、尺度和网格。
%   % 不是原始电影，也不是已去背景的电影。不要直接以原始 Y 替代 Ybg。
%
%   [dff,F0,fig] = plot_cnmfe_dff_svd(neuron,10,'F0',baseline);
%   % baseline: K x 1 或 K x T，必须与 neuron.C 的分量尺度一致。
%
% 可选参数:
%   'Signal'      'C' (默认) 或 'C_raw'
%   'Neurons'     [] 全部，或例如 1:35（矩阵行号）
%   'TimeRange'   [] 全部，或 [0 300] 秒
%   'CLim'        [] 自动，或 [0 12]
%   'Percentile'  50，背景基线分位数，范围 [0,100]
%   'WindowSec'   [] 全时段固定基线；正数为滑动窗口秒数
%   'OutputPrefix' '' 不保存；例如 'session01_dff' 保存 PNG 和 MAT
%
% 原理遵循 CNMF_E/@Sources2D/Sources2D.m 的
% extract_DF_F_endoscope: B_i=(A_i'*Ybg)/(A_i'*A_i), F0=quantile(B_i),
% dff=C/F0 或 C_raw/F0。C 已是去背景活动，不能再计算 (C-F0)/F0。
% 这是 CNMF-E 的背景参考 DF/F 估计，依赖背景估计质量；
% 仅有 C/C_raw 不能唯一恢复绝对荧光基线。不会以 eps 掩盖无效基线。
% 返回完整 K x T dff；Neurons/TimeRange 只影响显示，不影响计算。
% 不修改 neuron。使用 Ybg 分位数估计需要 prctile。

p=inputParser;
addRequired(p,'neuron',@(x)isstruct(x)||isobject(x));
addRequired(p,'Fs',@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'F0',[]); addParameter(p,'Ybg',[]);
addParameter(p,'Signal','C'); addParameter(p,'Neurons',[]);
addParameter(p,'TimeRange',[]); addParameter(p,'CLim',[]);
addParameter(p,'Percentile',50,@(x)isscalar(x)&&isfinite(x)&&x>=0&&x<=100);
addParameter(p,'WindowSec',[],@(x)isempty(x)||(isscalar(x)&&isfinite(x)&&x>0));
addParameter(p,'OutputPrefix','');
parse(p,neuron,Fs,varargin{:}); o=p.Results;
signal=validatestring(o.Signal,{'C','C_raw'});
C=double(full(neuron.(signal)));
assert(ismatrix(C)&&~isempty(C)&&isreal(C),'活动矩阵必须是非空实数 K x T 矩阵。');
[K,T]=size(C); assert(T>=2,'至少需要两帧。');
F0=[]; dff=[]; source='';
assert(isempty(o.F0)||isempty(o.Ybg),'F0 与 Ybg 只提供其中一个。');
if ~isempty(o.F0)
    F0=double(o.F0); source='user-supplied F0';
elseif ~isempty(o.Ybg)
    A=double(neuron.A); Ybg=o.Ybg;
    if ndims(Ybg)==3, Ybg=reshape(Ybg,[],size(Ybg,3)); end
    assert(size(A,2)==K && size(Ybg,1)==size(A,1) && size(Ybg,2)==T,...
        'Ybg 必须与 A 的像素数、C 的帧数一致。');
    normA=full(sum(A.^2,1))';
    assert(all(isfinite(normA)&normA>0),'A 存在空或无效空间分量。');
    % 分块投影，避免将完整背景电影复制为 double。
    B=zeros(K,T);
    for j=1:500:T
        jj=j:min(j+499,T);
        B(:,jj)=bsxfun(@rdivide,A'*double(Ybg(:,jj)),normA);
    end
    assert(all(isfinite(B(:))),'背景投影含 NaN/Inf，请先检查数据。');
    if isempty(o.WindowSec)
        F0=prctile(B,o.Percentile,2);
    else
        w=min(T,max(1,round(o.WindowSec*Fs)));
        F0=zeros(K,T); left=floor((w-1)/2); right=w-1-left;
        for j=1:T
            F0(:,j)=prctile(B(:,max(1,j-left):min(T,j+right)),o.Percentile,2);
        end
    end
    source='CNMF-E background projection percentile';
else
    if strcmp(signal,'C'), dfField='C_df'; else, dfField='C_raw_df'; end
    if hasvalue(neuron,dfField)
        dff=double(full(neuron.(dfField)));
        assert(isequal(size(dff),[K T]),'已有 DF/F 与活动矩阵尺寸不一致。');
        if hasvalue(neuron,'Df'), F0=double(neuron.Df); end
        source=['existing neuron.' dfField];
    elseif hasvalue(neuron,'Df')
        F0=double(neuron.Df); source='existing neuron.Df';
    else
        assert(isempty(o.WindowSec), ...
            'SVD 自动模式使用固定基线，请保持 WindowSec=[]。');
        F0=svd_baseline(neuron,o.Percentile);
        source='block SVD background: b*f+b0, projected percentile';
    end
end
if ~isempty(F0)
    if isvector(F0)&&numel(F0)==K, F0=F0(:); end
    assert(isequal(size(F0),[K 1])||isequal(size(F0),[K T]),'F0 必须是 K x 1 或 K x T。');
    bad = ~isfinite(F0) | F0 <= 0;

if any(bad(:))
    fprintf('\n基线来源：%s\n', source);

    % 每行对应一个神经元
    neuronID = (1:K)';
    nInvalid = sum(bad, 2);
    F0min = min(F0, [], 2);
    F0max = max(F0, [], 2);

    diagnostic = table(neuronID, F0min, F0max, nInvalid, ...
        'VariableNames', {'NeuronID','F0_min','F0_max','InvalidCount'});
    disp(diagnostic);

    warning('部分神经元的 F0 无效：对应 dff 将标为 NaN，并跳过绘图。');
end
end
% 保留全部神经元的原始行号，不修改 F0。
invalidNeuron = false(K,1);
if ~isempty(F0)
    invalidNeuron = any(~isfinite(F0) | F0 <= 0, 2);
end

if isempty(dff)
    % 无效基线不参与除法
    safeF0 = F0;
    safeF0(~isfinite(safeF0) | safeF0 <= 0) = NaN;
    dff = bsxfun(@rdivide, C, safeF0);
end

% 无效神经元整行标为 NaN
dff(invalidNeuron,:) = NaN;

% 其他神经元仍须通过有限值检查
validDff = dff(~invalidNeuron,:);
assert(all(isfinite(validDff(:))), ...
    '有效基线对应的活动仍含 NaN/Inf，请检查 C 或已有 DF/F。');

ids = o.Neurons;
if isempty(ids), ids = 1:K; end

assert(isvector(ids) && all(isfinite(ids) & ...
    ids >= 1 & ids <= K & ids == fix(ids)), ...
    'Neurons 必须是有效行号。');

ids = ids(~invalidNeuron(ids));
assert(~isempty(ids), '所选神经元均无有效基线，无法绘图。');
assert(isvector(ids)&&all(isfinite(ids)&ids>=1&ids<=K&ids==fix(ids)),'Neurons 必须是有效行号。');
t=(0:T-1)/Fs; keep=true(1,T);
if ~isempty(o.TimeRange)
    validateattributes(o.TimeRange,{'numeric'},{'real','finite','numel',2});
    assert(o.TimeRange(2)>o.TimeRange(1),'TimeRange 必须递增。');
    keep=t>=o.TimeRange(1)&t<=o.TimeRange(2);
end
assert(nnz(keep)>=2,'所选时间范围不足两帧。');
Z=dff(ids,keep); limits=o.CLim;
if isempty(limits)
    v=sort(Z(:)); lo=min(0,v(max(1,ceil(.01*numel(v)))));
    hi=v(max(1,ceil(.99*numel(v))));
    if hi<=lo, hi=lo+1; end
    limits=[lo hi];
end
validateattributes(limits,{'numeric'},{'real','finite','numel',2});
assert(limits(2)>limits(1),'CLim 必须递增。');
fig=figure('Color','w','Position',[150 250 1100 260]);
ax=axes('Parent',fig,'Position',[.065 .22 .84 .73]);
imagesc(ax,t(keep),1:numel(ids),Z); colormap(ax,parula(256)); caxis(ax,limits);
set(ax,'YDir','reverse','FontName','Arial','FontSize',11,...
    'TickDir','out','TickLength',[.004 .004],'Box','off','LineWidth',.8);
xlabel(ax,'Time (s)'); ylabel(ax,'Neuron ID');
step=max(1,ceil(numel(ids)/8)); ticks=unique([1 step:step:numel(ids)]);
set(ax,'YTick',ticks,'YTickLabel',arrayfun(@num2str,ids(ticks),'UniformOutput',false));
ylim(ax,[.5 numel(ids)+.5]);
cb=colorbar(ax); set(cb,'Units','normalized','Position',[.92 .22 .012 .73]);
cb.Label.String='\Delta F/F_0'; cb.Label.Interpreter='tex';
ax.Position=[.065 .22 .84 .73];
fprintf('DF/F source: %s\n',source);
if ~isempty(o.OutputPrefix)
    prefix=char(o.OutputPrefix);
    save([prefix '.mat'],'dff','F0','Fs','t','ids','source','-v7.3');
    print(fig,[prefix '.png'],'-dpng','-r300');
end
end

function tf=hasvalue(x,name)
if isstruct(x), tf=isfield(x,name); else, tf=isprop(x,name); end
if tf, tf=~isempty(x.(name)); end
end

function F0=svd_baseline(neuron,percentile)
% Equivalent to projecting the full SVD background, without allocating it.
assert(strcmpi(neuron.options.background_model,'svd'), ...
    '自动背景重建仅支持 background_model=svd。');
try
    patchPos=neuron.P.mat_data.patch_pos;
catch
    error(['无法读取分块位置。请确认 neuron.P.mat_file 指向的数据文件存在；' ...
        '若已移动，请先设置 neuron.P.mat_data=matfile(''数据文件的新路径'')。']);
end
A=double(neuron.A); [K,T]=size(neuron.C);
d1=neuron.options.d1; d2=neuron.options.d2;
assert(size(A,1)==d1*d2 && size(A,2)==K,'A 与 C/图像尺寸不一致。');
assert(iscell(patchPos)&&isequal(size(patchPos),size(neuron.b)), ...
    '分块位置与 b 的 cell 尺寸不一致。');
B=zeros(K,T); coverage=zeros(d1,d2); pixelIDs=reshape(1:d1*d2,d1,d2);
for k=1:numel(patchPos)
    pos=double(patchPos{k});
    assert(numel(pos)==4&&all(isfinite(pos))&&all(pos==fix(pos)), '无效分块坐标。');
    assert(pos(1)>=1&&pos(2)<=d1&&pos(3)>=1&&pos(4)<=d2 ...
        &&pos(2)>=pos(1)&&pos(4)>=pos(3),'分块坐标越界。');
    rr=pos(1):pos(2); cc=pos(3):pos(4);
    ix=pixelIDs(rr,cc); ix=ix(:); coverage(rr,cc)=coverage(rr,cc)+1;
    b=double(neuron.b{k}); f=double(neuron.f{k}); b0=double(neuron.b0{k});
    assert(size(b,1)==numel(ix)&&size(b,2)==size(f,1)&&size(f,2)==T ...
        &&numel(b0)==numel(ix),'第 %d 块背景尺寸不匹配。',k);
    Ap=A(ix,:);
    B=B+bsxfun(@plus,(Ap'*b)*f,Ap'*b0(:));
end
assert(all(coverage(:)==1),'分块必须无重叠且完整覆盖图像。');
normA=full(sum(A.^2,1))';
assert(all(isfinite(normA)&normA>0),'A 存在无效空间分量。');
B=bsxfun(@rdivide,B,normA);
assert(all(isfinite(B(:))),'背景投影含非有限值。');
if percentile==50, F0=median(B,2); else, F0=prctile(B,percentile,2); end
end
