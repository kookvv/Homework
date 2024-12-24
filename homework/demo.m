function varargout = demo(varargin)
% DEMO MATLAB code for demo.fig
%      DEMO, by itself, creates a new DEMO or raises the existing
%      singleton*.
%
%      H = DEMO returns the handle to a new DEMO or the handle to
%      the existing singleton*.
%
%      DEMO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DEMO.M with the given input arguments.
%
%      DEMO('Property','Value',...) creates a new DEMO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before demo_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to demo_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help demo

% Last Modified by GUIDE v2.5 24-Dec-2024 20:41:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @demo_OpeningFcn, ...
                   'gui_OutputFcn',  @demo_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before demo is made visible.
function demo_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to demo (see VARARGIN)

% Choose default command line output for demo
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes demo wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = demo_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
[filename,pathname]=uigetfile('*jpg','选择图片');
path=[pathname filename];
I=imread(path);
axes(handles.axes1)
imshow(I)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
% 如果是RGB图像，先转为灰度图像
if size(I, 3) == 3
    [height, width, ~] = size(I);
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 手动计算直方图
histogram = zeros(1, 256);  % 256个灰度级（0到255）
[height, width] = size(grayImage);
for i = 1:height
    for j = 1:width
        grayValue = round(grayImage(i, j));  % 获取每个像素的灰度值
        histogram(grayValue + 1) = histogram(grayValue + 1) + 1;  % 计数
    end
end

% 显示直方图
axes(handles.axes2);
bar(0:255, histogram, 'BarWidth', 1);
xlim([0 255]);
title('Gray Histogram');

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
% 如果是RGB图像，先转为灰度图像
if size(I, 3) == 3
    [height, width, ~] = size(I);
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 手动计算直方图
histogram = zeros(1, 256);
[height, width] = size(grayImage);
for i = 1:height
    for j = 1:width
        grayValue = round(grayImage(i, j));  % 获取每个像素的灰度值
        histogram(grayValue + 1) = histogram(grayValue + 1) + 1;  % 计数
    end
end

% 计算累计分布函数(CDF)
cdf = cumsum(histogram) / (height * width);  % 归一化CDF

% 直方图均衡化（映射灰度值）
equalizedImage = uint8(255 * cdf(round(grayImage) + 1));  % 映射到新的灰度值

% 在第二个坐标轴显示均衡化后的图像
axes(handles.axes2);
imshow(equalizedImage, []);

% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I  % I是目标图像

% 弹出文件选择对话框，选择参考图像
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp', 'Image Files (*.jpg, *.png, *.bmp)'}, 'Select Reference Image');
if isequal(filename, 0)
    disp('User selected Cancel');
    return;  % 如果用户取消了选择，退出函数
end

% 读取参考图像
targetImage = imread(fullfile(pathname, filename));  % 选择的参考图像
% 如果参考图像是RGB图像，先转为灰度图像
if size(targetImage, 3) == 3
    [height, width, ~] = size(targetImage);
    referenceImageGray = zeros(height, width);  % 初始化一个与参考图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = targetImage(i, j, 1);
            G = targetImage(i, j, 2);
            B = targetImage(i, j, 3);
            % 使用加权平均法转换为灰度图像
            referenceImageGray(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    referenceImageGray = targetImage;  % 如果已经是灰度图像，直接使用
end

% 如果目标图像是RGB图像，先转为灰度图像
if size(I, 3) == 3
    [height, width, ~] = size(I);
    targetImageGray = zeros(height, width);  % 初始化一个与目标图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            targetImageGray(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    targetImageGray = I;  % 如果已经是灰度图像，直接使用
end

% 计算目标图像和参考图像的灰度直方图
targetHist = zeros(1, 256);
referenceHist = zeros(1, 256);

[height, width] = size(targetImageGray);
for i = 1:height
    for j = 1:width
        targetValue = round(targetImageGray(i, j));  % 获取目标图像每个像素的灰度值
        referenceValue = round(referenceImageGray(i, j));  % 获取参考图像每个像素的灰度值
        targetHist(targetValue + 1) = targetHist(targetValue + 1) + 1;  % 计数目标图像的灰度值
        referenceHist(referenceValue + 1) = referenceHist(referenceValue + 1) + 1;  % 计数参考图像的灰度值
    end
end

% 计算目标图像和参考图像的累计分布函数 (CDF)
targetCDF = cumsum(targetHist) / (height * width);
referenceCDF = cumsum(referenceHist) / (height * width);

% 直方图匹配：为目标图像的每个像素值找到参考图像CDF最接近的值
mapping = zeros(1, 256);
for i = 0:255
    % 找到参考图像CDF中最接近目标图像CDF的灰度值
    [~, idx] = min(abs(referenceCDF - targetCDF(i + 1)));
    mapping(i + 1) = idx - 1;  % 映射关系
end

% 根据映射调整目标图像的像素值
matchedImage = uint8(zeros(height, width));
for i = 1:height
    for j = 1:width
        targetValue = round(targetImageGray(i, j));  % 获取目标图像每个像素的灰度值
        matchedImage(i, j) = mapping(targetValue + 1);  % 映射目标图像的灰度值
    end
end

% 显示直方图匹配后的图像
axes(handles.axes2);
imshow(matchedImage, []);

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I
% 检查图像是否是RGB图像（3通道）
if size(I, 3) == 3
    [height, width, ~] = size(I);
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 在第二个坐标轴显示灰度图像
axes(handles.axes2);
imshow(grayImage, []);


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I

% 检查图像是否是RGB图像（3通道）
if size(I, 3) == 3
    [height, width, ~] = size(I);
    % 将RGB图像转为灰度图像
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 对比度增强的线性变换
minI = min(grayImage(:));  % 图像中的最小灰度值
maxI = max(grayImage(:));  % 图像中的最大灰度值

% 设定线性变换的参数
alpha = 255 / (maxI - minI);  % 增强因子
beta = -alpha * minI;         % 偏移量

% 执行线性变换增强
enhancedImage = alpha * grayImage + beta;

% 防止出现溢出，确保图像值在[0, 255]范围内
enhancedImage(enhancedImage > 255) = 255;
enhancedImage(enhancedImage < 0) = 0;

% 将增强后的图像转换为 uint8 类型
enhancedImage = uint8(enhancedImage);

% 在第二个坐标轴显示增强后的图像
axes(handles.axes2);
imshow(enhancedImage, []);
   
% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I

% 获取输入图像的尺寸
[height, width, channels] = size(I);

% 设置缩放因子
scaleFactor = 0.8;

% 计算输出图像的尺寸
newHeight = round(height * scaleFactor);
newWidth = round(width * scaleFactor);

% 创建一张纯白色的图像，尺寸为原图的尺寸
whiteImage = 255 * ones(height, width, channels, 'uint8');

% 创建一个新的空图像，尺寸为新尺寸
resizedImage = zeros(newHeight, newWidth, channels, 'uint8');

% 计算双线性插值所需的坐标
for c = 1:channels  % 对于每个颜色通道（RGB）
    for i = 1:newHeight
        for j = 1:newWidth
            % 计算目标图像中的(i, j)对应原图像中的坐标
            origX = (j - 1) / scaleFactor + 1;
            origY = (i - 1) / scaleFactor + 1;
            
            % 计算原图像坐标的整数部分和小数部分
            x1 = floor(origX);
            y1 = floor(origY);
            x2 = min(x1 + 1, width);  % 防止越界
            y2 = min(y1 + 1, height); % 防止越界
            
            % 计算插值所需的权重
            dx = origX - x1;
            dy = origY - y1;
            
            % 获取四个相邻像素的值
            Q11 = double(I(y1, x1, c));  % 左上角
            Q12 = double(I(y2, x1, c));  % 左下角
            Q21 = double(I(y1, x2, c));  % 右上角
            Q22 = double(I(y2, x2, c));  % 右下角
            
            % 进行双线性插值
            resizedImage(i, j, c) = (1 - dx) * (1 - dy) * Q11 + dx * (1 - dy) * Q21 + (1 - dx) * dy * Q12 + dx * dy * Q22;
        end
    end
end

% 将缩放后的图像覆盖到白色背景图上
% 我们假设白色图像的尺寸足够容纳缩放后的图像
startX = round((width - newWidth) / 2);  % 计算缩放后的图像放置的位置（水平）
startY = round((height - newHeight) / 2);  % 计算缩放后的图像放置的位置（垂直）

% 将缩放后的图像粘贴到白色图像的中心位置
for c = 1:channels
    whiteImage(startY:startY+newHeight-1, startX:startX+newWidth-1, c) = resizedImage(:,:,c);
end

% 显示最终结果
axes(handles.axes2);
imshow(whiteImage);

% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I

% 获取原图像的尺寸
[height, width, channels] = size(I);

% 设置旋转角度（单位：度）
rotationAngle = 45;  % 可以根据需要修改角度

% 将角度转换为弧度
theta = deg2rad(rotationAngle);

% 计算旋转矩阵
rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% 计算旋转后图像的尺寸
newHeight = round(abs(height * cos(theta)) + abs(width * sin(theta)));
newWidth = round(abs(width * cos(theta)) + abs(height * sin(theta)));

% 创建一个全黑的图像，尺寸为旋转后的尺寸
rotatedImage = zeros(newHeight, newWidth, channels, 'uint8');

% 计算原图像的中心
centerX = round(width / 2);
centerY = round(height / 2);

% 计算旋转后图像的中心
newCenterX = round(newWidth / 2);
newCenterY = round(newHeight / 2);

% 遍历旋转后图像的每个像素
for i = 1:newHeight
    for j = 1:newWidth
        % 计算旋转后图像的像素位置相对于新图像中心的偏移
        offsetX = j - newCenterX;
        offsetY = i - newCenterY;
        
        % 逆向变换：将新图像的坐标映射到原图像坐标
        origCoords = rotationMatrix \ [offsetX; offsetY];
        
        origX = round(origCoords(1) + centerX);
        origY = round(origCoords(2) + centerY);
        
        % 如果映射到的原图像坐标在图像内，则赋值
        if origX >= 1 && origX <= width && origY >= 1 && origY <= height
            rotatedImage(i, j, :) = I(origY, origX, :);
        end
    end
end

% 显示旋转后的图像
axes(handles.axes2);
imshow(rotatedImage, []);

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I

% 检查图像是否是RGB图像（3通道）
if size(I, 3) == 3
    [height, width, ~] = size(I);
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 指数变换增强对比度
% 设置指数变换的常数 c 和指数 gamma
c = 1;    % 常数
gamma = 0.5;  % 指数，控制对比度增强的程度

% 对数变换后的图像初始化
expImage = zeros(size(grayImage));

% 计算指数变换
for i = 1:height
    for j = 1:width
        expImage(i, j) = c * (double(grayImage(i, j)) .^ gamma);  % 指数变换公式
    end
end

% 归一化指数变换后的图像到 [0, 255] 区间
expImage = uint8(255 * mat2gray(expImage));  % 使用 mat2gray 将结果归一化到 [0, 1] 然后乘以 255 得到 [0, 255] 的图像

% 在第二个坐标轴显示指数变换后的图像
axes(handles.axes2);
imshow(expImage, []);

% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global I

% 获取图像的大小
[height, width, channels] = size(I);

% 如果图像是RGB图像，先转换为灰度图像
if channels == 3
    grayImage = zeros(height, width);  % 初始化一个与输入图像同尺寸的灰度图像
    for i = 1:height
        for j = 1:width
            R = I(i, j, 1);
            G = I(i, j, 2);
            B = I(i, j, 3);
            % 使用加权平均法转换为灰度图像
            grayImage(i, j) = 0.2989 * R + 0.5870 * G + 0.1140 * B;
        end
    end
else
    grayImage = I;  % 如果已经是灰度图像，直接使用
end

% 对数变换增强对比度
% 设置对数变换的常数 C
C = 1;

% 计算图像的对数变换
logImage = zeros(size(grayImage));  % 初始化输出图像

% 对每个像素值进行对数变换，避免对数计算时出现负值，使用 max(grayImage(:)) 来确保像素值在对数变换时为正
maxPixelValue = max(grayImage(:));
for i = 1:height
    for j = 1:width
        logImage(i, j) = C * log(1 + double(grayImage(i, j)));  % 对数变换
    end
end

% 归一化对数变换后的图像到 [0, 255] 区间
logImage = uint8(255 * mat2gray(logImage));

% 在第二个坐标轴显示对数变换后的图像
axes(handles.axes2);
imshow(logImage, []);
