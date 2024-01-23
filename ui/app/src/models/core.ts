export interface File {
  path: string
  checksum: string
}

export interface Plate {
  barcode: string
  files: File[]
}

export interface Well {
  location: string
  plate: Plate
}

export interface Sample {
  well: Well
  pseudonymized_id: string
  plate: Plate
  sample_number?: string
  files: File[]
}

export interface Strain {
  name: string
}

export interface Panel {
  name: string
  strain: Strain
}

export interface Substrain {
  name: string
  strain: Strain
}

export interface SampleCount {
  substrain: Substrain
  sample: Sample
  plate: Plate
  aligned: number
  length: number
  rpkm: number
  rpkm_proportions: number
  outlier: boolean
  DP: number
  DP_threshold: number
  DP_status: string
  readnum_status: string
  readnum_threshold: number
  percentile_threshold: string
  normcounts?: number
  strain?: string
  pseudonymized_id?: string
}

export interface DataItem {
  panel: string
  value: number | string | null
  strain: string
}

export interface Metadata {
  plate: Plate
  sample: Sample
  well: Well
  prescriber: string
  order_date: string
  ent_date: string
  treatment_type: string
  data: Array<DataItem>
}
