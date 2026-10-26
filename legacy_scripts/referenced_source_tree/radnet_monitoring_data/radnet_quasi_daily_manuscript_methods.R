# ============================================================
# RadNet Gross Beta quasi-daily
# 按参考文章方法处理
#
# 核心原则：
#
# 1. 完全相同的：
#      loc_num
#      collect_start
#      collect_end
#    视为同一个滤膜采样周期；
#    若存在多条 Gross Beta，则取平均。
#
# 2. 每个不同滤膜周期独立保留。
#
# 3. 每个滤膜周期展开为：
#      collect_start <= date < collect_end
#
# 4. 不同滤膜周期即使日期发生重叠，
#    同一个 monitor + date 可以保留多行。
#
# 5. 绝不在 monitor-day 层面对不同滤膜周期
#    的 Gross Beta 再取平均。
# ============================================================

library(data.table)

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究/radnet 监测数据"

FILTER_FILE <- file.path(
  ROOT,
  "02_滤膜周期数据",
  "radnet_gross_beta_filter_period_2001_2024.csv.gz"
)

DAILY_FILE <- file.path(
  ROOT,
  "03_quasi-daily数据",
  "radnet_gross_beta_quasi_daily_2001_2024.csv.gz"
)

FINAL_FILE <- file.path(
  ROOT,
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)


# ============================================================
# 1. 读取滤膜周期数据
# ============================================================

x <- fread(FILTER_FILE)

x[
  ,
  collect_start := as.IDate(collect_start)
]

x[
  ,
  collect_end := as.IDate(collect_end)
]

x[
  ,
  result_date := as.IDate(result_date)
]


cat(
  "\n原始滤膜测量记录：",
  format(nrow(x), big.mark = ","),
  "\n",
  sep = ""
)

cat(
  "唯一 EPA loc_num 监测站：",
  uniqueN(x$loc_num),
  "\n",
  sep = ""
)

cat(
  "唯一 station + state：",
  uniqueN(x[, .(station, state_abbr)]),
  "\n",
  sep = ""
)


# ============================================================
# 2. 完全按照文章思想：
#
# 相同监测站 + 相同 collect_start + 相同 collect_end
# 才合并成同一个 filter period。
#
# 不同滤膜周期绝不合并。
# ============================================================

period <- x[
  !is.na(collect_start) &
    !is.na(collect_end) &
    !is.na(gross_beta_bq_m3),
  .(
    station =
      station[1],
    
    city_name =
      city_name[1],
    
    state_abbr =
      state_abbr[1],
    
    gross_beta_bq_m3 =
      mean(
        gross_beta_bq_m3,
        na.rm = TRUE
      ),
    
    csu_bq_m3 =
      if (
        all(is.na(csu_bq_m3))
      ) {
        NA_real_
      } else {
        mean(
          csu_bq_m3,
          na.rm = TRUE
        )
      },
    
    mdc_bq_m3 =
      if (
        all(is.na(mdc_bq_m3))
      ) {
        NA_real_
      } else {
        mean(
          mdc_bq_m3,
          na.rm = TRUE
        )
      },
    
    n_beta_measurements =
      .N,
    
    ana_num =
      paste(
        sort(
          unique(ana_num)
        ),
        collapse = ";"
      ),
    
    samp_num =
      paste(
        sort(
          unique(samp_num)
        ),
        collapse = ";"
      )
  ),
  by = .(
    loc_num,
    collect_start,
    collect_end
  )
]


# ============================================================
# 3. 计算真实采样天数
# ============================================================

period[
  ,
  collection_days :=
    as.integer(
      collect_end - collect_start
    )
]


# 生成一个唯一滤膜周期 ID
period[
  ,
  filter_period_id :=
    paste0(
      "L",
      loc_num,
      "_",
      format(
        collect_start,
        "%Y%m%d"
      ),
      "_",
      format(
        collect_end,
        "%Y%m%d"
      )
    )
]


setorder(
  period,
  loc_num,
  collect_start,
  collect_end
)


cat(
  "\n合并完全相同采样周期后：",
  format(
    nrow(period),
    big.mark = ","
  ),
  " 个 filter periods\n",
  sep = ""
)

cat(
  "其中同一周期包含 >1 条测量：",
  period[
    n_beta_measurements > 1,
    .N
  ],
  "\n",
  sep = ""
)

cat(
  "collection_days <= 0：",
  period[
    collection_days <= 0,
    .N
  ],
  "\n",
  sep = ""
)


# ============================================================
# 4. 只把有效 filter period 展开成 quasi-daily
#
# 文章规则：
#
# collect_start <= date < collect_end
#
# 即：
# collect_start 到 collect_end - 1
# ============================================================

period_valid <- period[
  collection_days > 0
]


quasi_daily <- period_valid[
  ,
  .(
    date =
      seq.Date(
        from =
          as.Date(
            collect_start[1]
          ),
        
        to =
          as.Date(
            collect_end[1]
          ) - 1,
        
        by = "day"
      )
  ),
  by = .(
    loc_num,
    station,
    city_name,
    state_abbr,
    
    filter_period_id,
    
    collect_start,
    collect_end,
    collection_days,
    
    gross_beta_bq_m3,
    csu_bq_m3,
    mdc_bq_m3,
    
    n_beta_measurements,
    ana_num,
    samp_num
  )
]


quasi_daily[
  ,
  date := as.IDate(date)
]


# ============================================================
# 5. 排序
#
# 注意：
# 不按 loc_num + date 聚合！
#
# 因此不同 filter period 即使覆盖同一天，
# 会自然保留成多行。
# ============================================================

setorder(
  quasi_daily,
  loc_num,
  date,
  collect_start,
  collect_end,
  filter_period_id
)


# ============================================================
# 6. 统计重叠 monitor-day
#
# 只打印，不删除，不平均。
# ============================================================

overlap <- quasi_daily[
  ,
  .(
    n_filter_periods =
      uniqueN(
        filter_period_id
      )
  ),
  by = .(
    loc_num,
    date
  )
][
  n_filter_periods > 1
]


cat(
  "\nquasi-daily period-day 总行数：",
  format(
    nrow(quasi_daily),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "唯一 monitor-day 数：",
  format(
    uniqueN(
      quasi_daily[
        ,
        .(
          loc_num,
          date
        )
      ]
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "存在多个 filter period 的 monitor-day：",
  format(
    nrow(overlap),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 7. 保存 quasi-daily
# ============================================================

fwrite(
  quasi_daily,
  DAILY_FILE,
  compress = "gzip",
  na = ""
)


# ============================================================
# 8. 正式模型表
#
# 文章式正式数据：
#
# observation unit =
# monitor × date × filter period
#
# 因此这里不再：
#
# group_by(loc_num, date)
# mean(gross_beta)
#
# 正式输出直接保留 filter_period_id。
# ============================================================

model_ready <- quasi_daily[
  ,
  .(
    loc_num,
    station,
    city_name,
    state_abbr,
    
    date,
    
    filter_period_id,
    
    collect_start,
    collect_end,
    collection_days,
    
    gross_beta_bq_m3,
    csu_bq_m3,
    mdc_bq_m3,
    
    n_beta_measurements,
    
    ana_num,
    samp_num
  )
]


setorder(
  model_ready,
  loc_num,
  date,
  collect_start,
  collect_end,
  filter_period_id
)


fwrite(
  model_ready,
  FINAL_FILE,
  compress = "gzip",
  na = ""
)


# ============================================================
# 9. 最终摘要
# ============================================================

cat(
  "\n",
  "============================================================\n",
  "RadNet quasi-daily 按文章方法处理完成\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nEPA loc_num 监测站数：",
  uniqueN(
    model_ready$loc_num
  ),
  "\n",
  sep = ""
)

cat(
  "filter periods：",
  format(
    nrow(period),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "quasi-daily period-day：",
  format(
    nrow(model_ready),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "唯一 monitor-day：",
  format(
    uniqueN(
      model_ready[
        ,
        .(
          loc_num,
          date
        )
      ]
    ),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "日期重叠 monitor-day：",
  format(
    nrow(overlap),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

cat(
  "\n重要：重叠 monitor-day 未平均、未删除，",
  "不同 filter period 分别保留。\n",
  sep = ""
)

cat(
  "\n正式输出：\n",
  FINAL_FILE,
  "\n",
  sep = ""
)

cat(
  "\n总体状态：完成\n"
)