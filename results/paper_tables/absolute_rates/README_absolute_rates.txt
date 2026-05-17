Absolute-rate supplementary tables
- included arm_type values: active
- events_total uses absolute_events when available, otherwise events
- AE-level table aggregates study × arm × molecule × window × AE first, deduplicating N per arm, then pools
- clinician_summary keeps molecule-level total participants visible and uses AE-specific denominators for each AE percentage
- clinician_wide reports each AE as events/N (%) for quick clinical reading
- global table uses 'any adverse event' rows when available
- if no any-AE rows are available, the global table reports denominator only; summing AE rows would double-count participants
