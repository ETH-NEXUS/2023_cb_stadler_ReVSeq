export interface TableDataSampleCounts {
  substrain: string
  strain: string
  aligned: number | string
  length: number | string
  rpkm: number | string
  rpkm_proportions: number | string
  normcounts?: number | string
}

export interface TableDataMetadata {
  pseudonymized_id: string
  well: string
  prescriber: string
  order_date: string | Date
  ent_date: string | Date
  treatment_type: string
}

export interface Column {
  name: string
  label: string
  field: string
  align: string
  sortable: boolean
}
