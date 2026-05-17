Absolute-rate supplementary tables
- included arm_type values: active
- events_total uses absolute_events when available, otherwise events
- AE-level table aggregates study × arm × molecule × window × AE first, deduplicating N per arm, then pools
- global table uses 'any adverse event' rows when available
- if no any-AE rows are available, the global table reports denominator only; summing AE rows would double-count participants
