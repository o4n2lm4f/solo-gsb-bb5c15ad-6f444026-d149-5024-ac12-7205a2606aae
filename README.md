# DNA 限制酶切模拟 API

基于 Python + FastAPI 的限制性内切酶消化模拟服务。输入一段仅含 A/C/G/T 的 DNA
序列（线性或环状）和若干酶（识别序列 + 切点偏移），返回按原链顺序排列的
片段序列、长度与边界。

## 启动（单条命令）

```bash
docker compose up --build
```

镜像构建时会自动安装全部依赖（`requirements.txt`），无需任何额外配置。

## 地址与端口

| 内容 | 地址 |
|---|---|
| API 基础地址 | http://localhost:8000 |
| 交互式文档 (Swagger UI) | http://localhost:8000/docs |
| 消化接口 | `POST http://localhost:8000/digest` |
| 健康检查 | `GET http://localhost:8000/health` |

默认端口 **8000**；如需修改，编辑 `docker-compose.yml` 中的端口映射
（如 `"8080:8000"` 将宿主端口改为 8080）。

## 请求格式

```json
{
  "sequence": "AAAAGAATTCTTTT",
  "topology": "linear",
  "enzymes": [
    { "name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1 }
  ]
}
```

| 字段 | 说明 |
|---|---|
| `sequence` | DNA 序列，仅允许 A/C/G/T（大小写均可，内部统一转为大写）；不能为空 |
| `topology` | `linear`（线性，默认）或 `circular`（环状） |
| `enzymes[].name` | 酶名称，非空 |
| `enzymes[].recognition` | 识别序列，仅允许 A/C/G/T，非空 |
| `enzymes[].cut_offset` | 切点偏移：从正向识别序列左端计，**含两端**，即 `0 <= cut_offset <= len(recognition)`。偏移 0 表示在识别序列左端之前切开，偏移等于识别长度表示在右端之后切开 |

## 切点规则

- 同时搜索识别序列及其**反向互补序列**；反向命中的切点按
  `识别长度 - 偏移` 换算（命中起点为 `i` 时切点为 `i + (k - offset)`）。
  回文识别序列（如 GAATTC）只搜索一次，同一位点不会重复切割。
- **环状**序列额外识别跨首尾的位点，切点对序列长度取模并去重。
- **线性**序列越界切点会被拒绝（不取模、不生效）。
- 切点 `p` 表示在第 `p-1` 与第 `p` 个碱基（0 起始）之间切开。
- 无切点时返回完整原序列（单片段）。

## 调用示例

```bash
curl -s -X POST http://localhost:8000/digest \
  -H "Content-Type: application/json" \
  -d '{
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}]
      }'
```

响应：

```json
{
  "topology": "linear",
  "sequence_length": 14,
  "cut_positions": [5],
  "fragment_count": 2,
  "fragments": [
    { "index": 0, "sequence": "AAAAG",     "length": 5, "start": 0, "end": 5  },
    { "index": 1, "sequence": "AATTCTTTT", "length": 9, "start": 5, "end": 14 }
  ],
  "enzymes": [
    { "name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1, "cuts": [5] }
  ]
}
```

片段按原链顺序排列；`start`/`end` 为原序列上的 0 起始坐标（`end` 不含）。
环状序列跨首尾的片段 `end <= start`，其序列为 `seq[start:] + seq[:end]`。

## 错误处理

非法输入返回 **422** 并说明原因，例如：

- 序列含非法碱基：`sequence contains invalid base(s): N, X (only A, C, G, T are allowed)`
- 空序列：`sequence must not be empty`
- 偏移越界：`cut_offset 7 is out of range for enzyme 'EcoRI': must be within [0, 6] ...`
- 非法拓扑：仅接受 `linear` / `circular`

## 本地开发与测试

```bash
pip install -r requirements.txt
uvicorn app.main:app --reload        # 开发模式启动
pip install httpx pytest
pytest tests/                        # 运行测试
```

## 项目结构

```
app/
  main.py      # FastAPI 入口与 /digest 端点
  digest.py    # 切点搜索与片段生成（纯逻辑，无框架依赖）
  schemas.py   # 请求/响应模型与校验
tests/test_digest.py
Dockerfile
docker-compose.yml
requirements.txt
```
