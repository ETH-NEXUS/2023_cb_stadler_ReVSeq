export interface Plate {
  barcode: string
}

export interface Well {
  location: string
  plate: Plate
}

export interface Sample {
  well: Well
  pseudoanonymized_id: string
  plate: Plate
  sample_number?: string
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
  qc_status: string
  coverage_threshold: number
  coverage: number
  normcounts?: number
  strain?: string
  pseudoanonymized_id?: string
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
